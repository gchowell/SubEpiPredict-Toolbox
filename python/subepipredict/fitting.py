"""Fitting layer: the multi-start fit of one candidate, the screen over the
number of sub-epidemics and the onset threshold, and the bootstrap.

Ports of ``fitSubepidemicModel.m``, ``fittingModifiedLogisticFunctionPatchABC.m``
and ``fittingModifiedLogisticFunctionPatchMultiple.m``.

The optimizer
-------------
MATLAB runs ``fmincon`` (algorithm ``sqp``) from every start point of a
``MultiStart``. The port runs SciPy's SLSQP, also a sequential quadratic
programming method, in the same coordinates (free parameters mapped to
[0, 1]; r, K and alpha on a log scale), on the same objective (log-likelihood
units divided by the number of observations), with the same finite-difference
step (1e-6) and function tolerance (1e-10), and from start points drawn the
same way. The two optimizers are different programs: they agree on a
well-determined minimum, and can end in different local minima of a
multimodal surface, as two MATLAB runs with different seeds do.
"""

from __future__ import annotations

import math
import os

import numpy as np
from scipy.optimize import Bounds, minimize
from scipy.special import gammaln

from .aicc import get_nparams, getAICc
from .errors import SubEpiPredictError, subepidemicWarning
from .identifiability import subepidemicBoundSummary, subepidemicIdentifiability, subepidemicPrintPII
from .likelihood import subepidemicObjective
from .model import simulateSubepidemic, subepidemicOdeOptions, unpackSubepidemicParameters
from .noise import AddPoissonError, getMeanVarianceRatio
from .search import (
    _map_toU,
    _map_toZ,
    subepidemicOrderLabels,
    subepidemicSearchSpace,
    subepidemicStartPoints,
    subepidemicUnitMap,
    unique_rows_stable,
)
from .settings import resolve_fit_settings, subepidemicFitSeries
from .utils import Struct, get_rng, mat2str, movmean, num2str

_REALMIN = np.finfo(float).tiny
_PENALTY = 1e100  # objective of an invalid parameter vector, log-likelihood units

# fmincon option names accepted in settings.optimOptions.
_DEFAULT_OPTIM = {
    "Algorithm": "sqp",
    "StepTolerance": 1e-10,
    "FunctionTolerance": 1e-10,
    "OptimalityTolerance": 1e-8,
    "FiniteDifferenceType": "forward",
    "FiniteDifferenceStepSize": 1e-6,
    "MaxFunctionEvaluations": 20000,
    "MaxIterations": 20000,
}

_CONTEXT_FIELDS = ("timevect", "ydata", "I0", "npatches", "onset_fixed", "onset_thr", "flag1", "method1", "odeOptions")


# <============================================================================>
# <============================ Scaled objective ==============================>
# <============================================================================>


class _Objective:
    """The objective the optimizer sees, as a picklable callable.

    ``natural(u)`` is the objective in log-likelihood units at the point ``u``
    of the unit box (1e100 for an invalid parameter vector); calling the
    object returns it divided by the number of observations.
    """

    def __init__(self, context, map_, nObs, dataConstant):
        self.context = context
        self.map = map_
        self.nObs = nObs
        self.dataConstant = dataConstant

    def natural(self, u):
        v, detail = subepidemicObjective(_map_toZ(u, self.map), self.context)
        if not detail.success:
            return _PENALTY
        method1 = self.context["method1"]
        if method1 == 0:
            value = 0.5 * self.nObs * math.log(max(v, _REALMIN))
        elif method1 == 6:
            value = self.nObs * math.log(max(v, _REALMIN))
        else:
            value = v + self.dataConstant
        if not math.isfinite(value):
            return _PENALTY
        return value

    def __call__(self, u):
        return self.natural(u) / self.nObs


class _EvaluationLimit(Exception):
    pass


def _local_solve(fs, u0, nFree, optim, penalty):
    """One bounded local minimization in the unit box (stand-in for ``fmincon``/sqp).

    Returns ``(u, f, exitflag)`` with MATLAB-like exit flags: 1 converged,
    0 ran out of iterations or evaluations, -2 no usable point. ``penalty`` is
    the value ``fs`` returns for an invalid parameter vector.

    The gradient is taken by finite differences with the step of fmincon. A
    difference is not taken across the edge of the valid region: when the
    forward point is invalid the backward point is used.

    SLSQP is run in rounds, each from the end point of the previous one, until
    the objective no longer decreases. SLSQP can stop early on a badly scaled
    stretch of the surface; the rounds are what the step-size test of fmincon
    achieves.
    """
    h = float(optim["FiniteDifferenceStepSize"])
    central = str(optim.get("FiniteDifferenceType", "forward")).lower() == "central"
    maxfev = int(optim["MaxFunctionEvaluations"])
    state = {"n": 0, "best_f": math.inf, "best_u": None, "base_u": None, "base_f": None}
    scale = [1.0]

    def F(u):
        if state["n"] >= maxfev:
            raise _EvaluationLimit()
        state["n"] += 1
        v = fs(u)
        if v < state["best_f"]:
            state["best_f"] = v
            state["best_u"] = np.array(u, dtype=float, copy=True)
        return v

    def base(u):
        if state["base_u"] is None or not np.array_equal(u, state["base_u"]):
            state["base_f"] = F(u)
            state["base_u"] = np.array(u, dtype=float, copy=True)
        return state["base_f"]

    def gradient(u):
        f0 = base(u)
        g = np.zeros(nFree)
        for j in range(nFree):
            canF = u[j] + h <= 1.0
            canB = u[j] - h >= 0.0
            fp = fm = None
            if canF:
                up = u.copy()
                up[j] += h
                fp = F(up)
                if fp >= penalty:
                    fp = None
            if canB and (fp is None or central):
                um = u.copy()
                um[j] -= h
                fm = F(um)
                if fm >= penalty:
                    fm = None
            if f0 >= penalty:
                continue
            if fp is not None and fm is not None:
                g[j] = (fp - fm) / (2 * h)
            elif fp is not None:
                g[j] = (fp - f0) / h
            elif fm is not None:
                g[j] = (f0 - fm) / h
        return g

    def fun(u):
        return base(np.asarray(u, dtype=float)) / scale[0]

    def jac(u):
        return gradient(np.asarray(u, dtype=float)) / scale[0]

    bounds = Bounds(np.zeros(nFree), np.ones(nFree), keep_feasible=True)
    max_rounds = 5
    iters_left = int(optim["MaxIterations"])
    options = {"ftol": float(optim["FunctionTolerance"]), "disp": False}

    u = np.minimum(np.maximum(np.asarray(u0, dtype=float).ravel(), 0.0), 1.0)
    exitflag = -2
    try:
        f_prev = base(u)
        if f_prev >= penalty:
            # An invalid start point. Look at its neighbours (the points a
            # finite-difference gradient visits): continue from the best valid
            # one, or report the flat penalty plateau as fmincon does, as a
            # converged run whose value is the penalty.
            gradient(u)
            if state["best_f"] < penalty:
                u = state["best_u"].copy()
                f_prev = base(u)
            else:
                return u, float(f_prev), 1
        for _ in range(max_rounds):
            if iters_left <= 0:
                break
            scale[0] = max(1.0, abs(f_prev))
            options["maxiter"] = int(min(iters_left, 200))
            res = minimize(fun, u, jac=jac, method="SLSQP", bounds=bounds, options=options)
            iters_left -= int(getattr(res, "nit", 0))
            if state["best_u"] is None:
                break
            u_new = state["best_u"].copy()  # best point visited in this run so far
            f_new = float(state["best_f"])
            if res.status == 9:
                exitflag = 0  # iteration limit
            elif res.success or res.status == 8:
                # status 8 ('positive directional derivative for linesearch') is how
                # SLSQP reports that no further descent is possible at the precision
                # of the finite-difference gradient.
                exitflag = 1
            else:
                exitflag = -2
            improved = f_prev - f_new > 1e-9 * max(1.0, abs(f_new))
            if f_new <= f_prev:
                u = u_new
                f_prev = f_new
            if exitflag < 0 or not improved:
                if exitflag == 0 and not improved:
                    exitflag = 1  # stalled at the precision of the objective
                break
    except _EvaluationLimit:
        exitflag = 0

    if state["best_u"] is None or not math.isfinite(state["best_f"]):
        return u, math.inf, -2
    return np.minimum(np.maximum(state["best_u"], 0.0), 1.0), float(state["best_f"]), exitflag


def _solve_task(task):
    fs, u0, nFree, optim, penalty = task
    try:
        return _local_solve(fs, u0, nFree, optim, penalty)
    except Exception:  # a failure inside the local solver: this start gives no solution
        return None


# ---- process pool for settings.useParallel -----------------------------------

_POOL = None


def _get_pool():
    """The process pool of the toolbox, started on first use (like ``parpool``)."""
    global _POOL
    if _POOL is None:
        import atexit
        from concurrent.futures import ProcessPoolExecutor

        _POOL = ProcessPoolExecutor()
        atexit.register(shutdown_pool)
    return _POOL


def shutdown_pool():
    """Stop the process pool started by ``useParallel = True`` (it is also stopped at exit)."""
    global _POOL
    if _POOL is not None:
        try:
            _POOL.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
        _POOL = None


def _multistart(fs, U0, nFree, optim, penalty, pool=None):
    """Run the local solver from every start point (stand-in for ``MultiStart``).

    Returns ``(x, fval, exitflag, minima)``. ``minima`` holds the DISTINCT
    local minima of the runs with a positive exit flag, best first, as the
    solutions output of MultiStart does (two runs are one minimum when their
    objective values agree within 1e-6 relative and their points within 1e-6).
    The overall exit flag follows MultiStart: 1 all runs converged, 2 some
    did, 0 none converged but a run stopped at its iteration limit, -8 no
    solution.
    """
    tasks = [(fs, U0[k], nFree, optim, penalty) for k in range(U0.shape[0])]
    if pool is not None and len(tasks) > 1:
        results = list(pool.map(_solve_task, tasks))
    else:
        results = [_solve_task(task) for task in tasks]
    runs = [Struct(X=r[0], Fval=r[1], Exitflag=r[2]) for r in results if r is not None]

    if not runs:
        return None, math.inf, -10, []

    converged = sorted((r for r in runs if r.Exitflag > 0), key=lambda r: r.Fval)
    minima = []
    for r in converged:
        for q in minima:
            if abs(q.Fval - r.Fval) <= 1e-6 * max(1.0, abs(q.Fval)) and np.linalg.norm(q.X - r.X) <= 1e-6 * max(
                1.0, np.linalg.norm(q.X)
            ):
                break
        else:
            minima.append(r)

    if converged:
        best = converged[0]
        flag = 1 if len(converged) == len(runs) else 2
        return best.X, best.Fval, flag, minima

    stalled = [r for r in runs if r.Exitflag == 0]
    if stalled:
        best = min(stalled, key=lambda r: r.Fval)
        return best.X, best.Fval, 0, minima

    return None, math.inf, -8, minima


def _localGradient(natU, uP, h):
    """Projected central-difference gradient in unit coordinates, log-likelihood units.

    One-sided at a bound. Components that push into an active bound are set
    to zero. Returns ``(g, gInf)``.
    """
    n = np.size(uP)
    g = np.full(n, np.nan)
    gInf = math.nan
    try:
        uP = np.asarray(uP, dtype=float).ravel()
        g = np.zeros(n)
        fP = None
        for j in range(n):
            e = np.zeros(n)
            e[j] = h
            if uP[j] + h <= 1 and uP[j] - h >= 0:
                g[j] = (natU(uP + e) - natU(uP - e)) / (2 * h)
            else:
                if fP is None:
                    fP = natU(uP)
                if uP[j] + h <= 1:
                    g[j] = (natU(uP + e) - fP) / h
                else:
                    g[j] = (fP - natU(uP - e)) / h
        g[(uP <= 1e-8) & (g > 0)] = 0  # at the lower bound, descent would leave the box
        g[(uP >= 1 - 1e-8) & (g < 0)] = 0  # at the upper bound, likewise
        gInf = float(np.max(np.abs(g))) if n else 0.0
    except Exception:
        pass  # diagnostics only
    return g, gInf


# <============================================================================>
# <============================ fitSubepidemicModel ===========================>
# <============================================================================>


def fitSubepidemicModel(context, z0, LB, UB, isLog, starts=None, numstartpoints=0, settings=None, rng=None):
    """Multi-start fit of one candidate in a scaled search space.

    Port of ``fitSubepidemicModel.m``.

    Parameters
    ----------
    context : dict / Struct
        Fit context of :func:`subepidemicObjective`: ``timevect``, ``ydata``,
        ``I0``, ``npatches``, ``onset_fixed``, ``onset_thr``, ``flag1``,
        ``method1``, ``odeOptions``.
    z0, LB, UB : seed and bounds, layout ``[r, p, a, K, alpha, d]`` with
        ``4*n + 2`` entries (:func:`subepidemicSearchSpace`). A parameter with
        ``LB == UB`` is fixed at that value.
    isLog : same length; coordinates searched on a log scale
    starts : None or a matrix of start points, one per row in original units,
        used for the FIRST multi-start pass. With None the function draws
        ``numstartpoints`` Latin-hypercube points, jittered copies of ``z0``
        and ``z0`` itself.
    numstartpoints : Latin-hypercube points per generated start set
    settings : :func:`subepidemicFitSettings` (optional)
    rng : random stream for generated start points (default: toolbox stream)

    Returns
    -------
    ``(P, fval, fitInfo)``
        ``P`` the fitted vector in original units (NaN when the fit failed);
        ``fval`` the objective at ``P`` in the legacy units (SSE, negative
        log-likelihood without its data constant, or SAD), the number
        :func:`getAICc` expects (``inf`` when the fit failed); ``fitInfo``
        diagnostics. ``fitInfo.success`` must be checked.

    What the optimizer sees: only the free parameters, each mapped to [0, 1];
    the objective in log-likelihood units divided by the number of
    observations (least squares ``(n/2) log(SSE)``, absolute deviations
    ``n log(SAD)``, count likelihoods the legacy objective plus
    ``sum(gammaln(floor(y)+1))``). An invalid parameter vector gets a penalty
    of 1e100; every local minimum is evaluated again in original units and
    must be valid before it is accepted.

    Diagnostics in ``fitInfo``: ``success``, ``message``, ``exitflag``,
    ``attempts``, ``nActive`` (sub-epidemics active by the end of the
    calibration grid at ``P``), ``startAgreement`` (fraction of converged
    local minima within ``settings.agreementTol`` of the best one),
    ``nLocalRuns``, ``nConvergedRuns``, ``gradient`` and ``gradientInfNorm``
    (projected central-difference gradient at ``P``, log-likelihood units per
    unit coordinate), ``polish`` (Struct; the Nelder-Mead polish runs only
    when ``gradientInfNorm`` exceeds ``settings.polishGradientTol``),
    ``atLowerBound``, ``atUpperBound`` (ZERO-BASED indices into ``P``),
    ``fallback`` (when ``success`` is False: the best valid point any local
    run ended at, or None), ``freeIndex`` (zero-based), ``nFree``, ``isLog``,
    ``optimizerOptions``, ``useParallel``, ``objectiveUnits``.
    """
    settings = resolve_fit_settings(settings)
    rng = get_rng(rng)
    numstartpoints = 0 if numstartpoints is None else numstartpoints
    if not (np.ndim(numstartpoints) == 0 and math.isfinite(numstartpoints) and numstartpoints >= 0
            and numstartpoints == int(numstartpoints)):  # fmt: skip
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidStartCount", "numstartpoints must be a nonnegative integer."
        )
    numstartpoints = int(numstartpoints)
    for name in _CONTEXT_FIELDS:
        if name not in context:
            raise SubEpiPredictError("SubEpiPredict:FitContext", "The fit context has no field '%s'." % name)
    if context["method1"] not in (0, 1, 3, 4, 5, 6):
        raise SubEpiPredictError("SubEpiPredict:InvalidMethod", "method1 must be 0, 1, 3, 4, 5 or 6.")
    nz = 4 * int(context["npatches"]) + 2
    if any(np.size(v) != nz for v in (z0, LB, UB, isLog)):
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Seed, bounds and scale flags must all have 4*npatches+2 = %d entries." % nz,
        )
    map_ = subepidemicUnitMap("make", dict(z0=z0, LB=LB, UB=UB, isLog=isLog))
    z0 = map_.z0
    nFree = map_.nFree
    ydata = np.asarray(context["ydata"], dtype=float).ravel()
    nObs = ydata.size
    if nObs < 2 or np.any(~np.isfinite(ydata)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidData", "The fitted series must hold at least two finite real values."
        )

    # Objective in log-likelihood units. The data constant restores the term the
    # count objectives drop, so their value is a full negative log-likelihood.
    dataConstant = 0.0
    if context["method1"] in (1, 3, 4, 5):
        dataConstant = float(np.sum(gammaln(np.floor(ydata) + 1)))
    fs = _Objective(context, map_, nObs, dataConstant)  # what the optimizer minimizes
    natU = fs.natural
    penalty = 0.5 * _PENALTY / nObs

    P = np.full(nz, np.nan)
    fval = math.inf
    fitInfo = Struct(
        success=False, message="", exitflag=math.nan, attempts=0, attemptExitflags=[], nActive=0,
        freeIndex=map_.free, nFree=nFree, isLog=map_.isLog,
        objectiveUnits="log-likelihood units; the optimizer minimizes them divided by the number of observations",
        useParallel=False, nLocalRuns=0, nConvergedRuns=0, startAgreement=math.nan,
        gradient=np.zeros(0), gradientInfNorm=math.nan,
        atLowerBound=np.zeros(0, dtype=int), atUpperBound=np.zeros(0, dtype=int), fallback=None,
    )  # fmt: skip
    fitInfo.attemptMessages = []
    fitInfo.optimizerOptions = Struct()
    fitInfo.polish = Struct(
        enabled=settings.polish, run=False, improved=False, exitflag=math.nan, confirmExitflag=math.nan,
        objectiveBefore=math.nan, objectiveAfter=math.nan, gradientBefore=math.nan,
        gradientTol=settings.polishGradientTol, budget=settings.polishBudget,
    )  # fmt: skip

    accepted = False
    chosenDetail = None
    fallback = None

    if nFree == 0:
        # nothing to estimate: evaluate the fixed vector
        value, detail = subepidemicObjective(z0, context)
        if detail.success:
            P = z0.copy()
            fval = value
            chosenDetail = detail
            accepted = True
            fitInfo.exitflag = 1
        else:
            fitInfo.message = detail.message
    else:
        # ---- optimizer settings (defaults, then the user's overrides) --------
        optim = dict(_DEFAULT_OPTIM)
        for key, value in settings.optimOptions.items():
            if key not in optim and key not in ("Display",):
                raise SubEpiPredictError(
                    "SubEpiPredict:FitSettings",
                    "optimOptions has the unknown optimizer option '%s'. The options are %s."
                    % (key, ", ".join(sorted(_DEFAULT_OPTIM))),
                )
            optim[key] = value
        fdStep = float(np.ravel(optim["FiniteDifferenceStepSize"])[0])
        optim["FiniteDifferenceStepSize"] = fdStep
        fitInfo.optimizerOptions = Struct({k: optim[k] for k in _DEFAULT_OPTIM})

        # ---- multi-start, serial or parallel ---------------------------------
        pool = None
        if settings.useParallel:
            try:
                pool = _get_pool()
                fitInfo.useParallel = True
            except Exception as ME:
                subepidemicWarning(
                    "SubEpiPredict:ParallelUnavailable",
                    "Parallel pool unavailable (%s); running the multi-start serially.", ME,
                )  # fmt: skip
                pool = None

        nJitter = min(3, max(1, numstartpoints))
        for attempt in range(1, int(settings.fitMaxAttempts) + 1):
            if attempt == 1 and starts is not None and np.size(starts) > 0:
                startSet = np.asarray(starts, dtype=float)
                if startSet.ndim == 1:
                    startSet = startSet.reshape(1, -1)
                if startSet.shape[1] != nz or np.any(~np.isfinite(startSet)):
                    raise SubEpiPredictError(
                        "SubEpiPredict:ParameterDimensionMismatch",
                        "Start points must be finite, one per row, with %d columns." % nz,
                    )
            else:
                startSet = subepidemicStartPoints(map_, numstartpoints, nJitter, z0, rng)
            U0, _ = unique_rows_stable(_map_toU(startSet, map_))
            fitInfo.attempts = attempt
            try:
                candidate, value, flag, minima = _multistart(fs, U0, nFree, optim, penalty, pool)
                message = ""
            except Exception as ME:
                if pool is not None:
                    # a broken pool: run this pass, and the rest of the fit, serially
                    subepidemicWarning(
                        "SubEpiPredict:ParallelUnavailable",
                        "Parallel pool unavailable (%s); running the multi-start serially.", ME,
                    )  # fmt: skip
                    shutdown_pool()
                    pool = None
                    fitInfo.useParallel = False
                    candidate, value, flag, minima = _multistart(fs, U0, nFree, optim, penalty, None)
                    message = ""
                else:
                    flag, candidate, value, minima = -10, None, math.inf, []
                    message = "%s: %s" % (type(ME).__name__, ME)
            fitInfo.attemptExitflags.append(flag)
            fitInfo.attemptMessages.append(message)
            fitInfo.exitflag = flag
            fitInfo.nLocalRuns = U0.shape[0]
            fitInfo.nConvergedRuns = 0

            # Start agreement: how many converged local minima reproduce the best one.
            if minima:
                fvs = np.array([m.Fval for m in minima]) * nObs  # log-likelihood units
                okRun = np.isfinite(fvs)
                fitInfo.nConvergedRuns = int(np.sum(okRun))
                if np.any(okRun):
                    bestOK = float(np.min(fvs[okRun]))
                    fitInfo.startAgreement = float(
                        np.sum(okRun & (fvs <= bestOK + settings.agreementTol)) / np.sum(okRun)
                    )

            # A penalty plateau can appear to converge. Evaluate each local minimum
            # again in original units and accept the best one that is valid.
            candidates, values, flags = [], [], []
            if flag > 0 and candidate is not None and math.isfinite(value):
                candidates.append(_map_toZ(candidate, map_))
                values.append(value)
                flags.append(flag)
            for m in minima:
                if m.Exitflag > 0 and math.isfinite(m.Fval):
                    candidates.append(_map_toZ(m.X, map_))
                    values.append(m.Fval)
                    flags.append(m.Exitflag)
            for j in np.argsort(values, kind="stable"):
                trial = candidates[j]
                value, detail = subepidemicObjective(trial, context)
                if detail.success:
                    P = trial.copy()
                    fval = value
                    chosenDetail = detail
                    fitInfo.exitflag = flags[j]
                    accepted = True
                    break
                fitInfo.message = detail.message
            if accepted:
                break

            # Not accepted: remember the best valid end point for the fallback record.
            if candidate is not None and np.size(candidate) == nFree and np.all(np.isfinite(candidate)):
                trial = _map_toZ(candidate, map_)
                value, detail = subepidemicObjective(trial, context)
                if detail.success and (fallback is None or value < fallback.fval):
                    fallback = Struct(P=trial, fval=value, nActive=detail.nActive)
            # An unsuccessful pass is retried with fresh start points.

        # ---- stationarity check, then (only if needed) the polish ------------
        if accepted:
            g, gInf = _localGradient(natU, _map_toU(P, map_)[0], fdStep)
            fitInfo.gradient = g
            fitInfo.gradientInfNorm = gInf
            fitInfo.polish.gradientBefore = gInf
        if accepted and settings.polish and not (
            math.isfinite(fitInfo.gradientInfNorm) and fitInfo.gradientInfNorm <= settings.polishGradientTol
        ):
            fitInfo.polish.run = True
            fitInfo.polish.objectiveBefore = fval
            uP = _map_toU(P, map_)[0]
            budget = int(math.ceil(settings.polishBudget * nFree))

            def boxPenalty(u):
                # Nelder-Mead objective: the scaled objective inside the unit box, a penalty outside.
                u = np.asarray(u, dtype=float).ravel()
                if np.any(u < 0) or np.any(u > 1):
                    return _PENALTY
                return fs(u)

            try:
                res = minimize(
                    boxPenalty, uP, method="Nelder-Mead",
                    options={"xatol": 1e-6, "fatol": 1e-6, "maxfev": budget, "maxiter": budget, "disp": False},
                )  # fmt: skip
                uRaw = np.asarray(res.x, dtype=float).ravel()
                polishFlag = 1 if res.success else 0
                fitInfo.polish.exitflag = polishFlag
                uTry = np.minimum(np.maximum(uRaw, 0.0), 1.0)
                zTry = _map_toZ(uTry, map_)
                vTry, dTry = subepidemicObjective(zTry, context)  # validity check in original units
                if dTry.success and vTry < fval:
                    P = zTry
                    fval = vTry
                    chosenDetail = dTry
                    fitInfo.polish.improved = True
                    # A polish that improved but did not converge (budget used up, or
                    # an optimum on a bound) is handed to a bounded sqp run from the
                    # polished point; its exit flag describes the returned vector.
                    if polishFlag != 1 or np.any(uRaw < 0) or np.any(uRaw > 1):
                        try:
                            uC, _, flagC = _local_solve(fs, uTry, nFree, optim, penalty)
                            fitInfo.polish.confirmExitflag = flagC
                            zC = _map_toZ(uC, map_)
                            vC2, dC = subepidemicObjective(zC, context)
                            if flagC > 0 and dC.success and vC2 <= fval:
                                P = zC
                                fval = vC2
                                chosenDetail = dC
                                fitInfo.exitflag = flagC
                        except Exception:
                            pass  # keep the polished point
            except Exception:
                pass  # keep the sqp point if the polish fails for any reason
            fitInfo.polish.objectiveAfter = fval
            if fitInfo.polish.improved:  # check stationarity again at the returned point
                g, gInf = _localGradient(natU, _map_toU(P, map_)[0], fdStep)
                fitInfo.gradient = g
                fitInfo.gradientInfNorm = gInf

    if accepted:
        fitInfo.success = True
        fitInfo.nActive = chosenDetail.nActive
        fitInfo.message = "Accepted a converged local minimum with a valid trajectory."
        if nFree > 0:
            uP = _map_toU(P, map_)[0]
            fitInfo.atLowerBound = map_.free[uP <= 1e-6]
            fitInfo.atUpperBound = map_.free[uP >= 1 - 1e-6]
    else:
        # The seed is the last resort of the fallback record, never an accepted fit.
        value, detail = subepidemicObjective(z0, context)
        if detail.success and (fallback is None or value < fallback.fval):
            fallback = Struct(P=z0.copy(), fval=value, nActive=detail.nActive)
        fitInfo.fallback = fallback
        if not fitInfo.message:
            fitInfo.message = "No converged, valid local minimum. Seed diagnostic: " + detail.message
        P = np.full(nz, np.nan)
        fval = math.inf
    return P, fval, fitInfo


# <============================================================================>
# <============== Screen over sub-epidemic numbers and thresholds =============>
# <============================================================================>

THRESHOLD_DIAGNOSTICS_COLUMNS = ["nCandidate", "onset_thr", "nActive", "retained", "fval"]
FIT_DIAGNOSTICS_COLUMNS = ["nCandidate", "onset_thr", "success", "exitflag", "startAgreement", "gradientInfNorm",
                           "polishRun", "polishImproved", "nLocalRuns", "nConvergedRuns", "nAtBound"]  # fmt: skip


def onsetThresholdGrid(incidence, smoothfactor1):
    """Grid of onset thresholds C_thr of the screen.

    Equal-width discretization of the cumulative curve of the moving average
    of the calibration incidence: 0 followed by ``n - 1`` equally spaced
    values from its first value towards its last (the last is left out), for
    ``n`` data points.
    """
    incidence = np.asarray(incidence, dtype=float).ravel()
    cumcurve1 = np.cumsum(movmean(incidence, smoothfactor1))
    onset_thrs = np.linspace(cumcurve1[0], cumcurve1[-1], incidence.size)
    return np.concatenate(([0.0], onset_thrs[:-1]))


def _prepare_window(data1, epidemic_period):
    """Rows ``epidemic_period`` (one-based) of ``data1``; a leading zero is dropped."""
    data1 = np.asarray(data1, dtype=float)
    rows = np.asarray(epidemic_period, dtype=int).ravel() - 1
    data1 = data1[rows, :]
    if data1[0, 1] == 0:
        data1 = data1[1:, :]
    return data1


def fittingModifiedLogisticFunctionPatchABC(
    datafilename1, data1, DT, epidemic_period, M, flagX, numstartpoints, *,
    method1, dist1, npatches_fixed, onset_fixed, smoothfactor1, calibrationperiod1,
    fit_settings=None, output_dir="./output", save=True, plot=True, verbose=True, rng=None,
):  # fmt: skip
    """Rank n-sub-epidemic models over the number of sub-epidemics and the onset threshold.

    Port of ``fittingModifiedLogisticFunctionPatchABC.m``. For every number of
    sub-epidemics from 1 to ``npatches_fixed`` and every onset threshold of
    the grid (:func:`onsetThresholdGrid`), the model is fitted by
    :func:`fitSubepidemicModel` and its AICc recorded. A candidate is kept
    only when all its sub-epidemics are active by the end of the calibration
    period. The candidates are sorted by AICc.

    Parameters
    ----------
    datafilename1 : name of the adjusted data file of the area (used in the
        name of the saved file)
    data1 : n-by-2 array, time index and incidence
    DT : time step (not used by the fit; kept as in MATLAB)
    epidemic_period : ONE-BASED rows of ``data1`` to use
    M : number of bootstrap realizations (not used here; kept as in MATLAB)
    flagX : growth model, 0 to 5
    numstartpoints : Latin-hypercube start points per fit
    method1, dist1, npatches_fixed, onset_fixed, smoothfactor1, calibrationperiod1 :
        the settings MATLAB reads from global variables
    fit_settings : fit settings (see :func:`subepidemicFitSettings`)
    output_dir, save : folder of the saved ``ABC-ensem-...mat`` file
    plot : draw figure 100 (best fit of the screen)
    rng : random stream of the start points

    Returns
    -------
    ``(RMSES, PS, npatches, onset_thr, P, screen)``
        ``RMSES``: one row per retained candidate, ``[number of
        sub-epidemics, onset threshold, AICc]``, best first. ``PS``: the
        fitted vector of each candidate in the first ``4*n + 2`` columns of
        its row (the rest are zeros). ``npatches``, ``onset_thr``, ``P``: the
        best candidate. ``screen``: a Struct with everything else the MATLAB
        function saves (``thresholdDiagnostics``, ``fitDiagnostics``,
        ``candidateFitInfo``, ``relativelik_i``, ``timevect``, ``I0``,
        ``data``, ``ydata``, ``bestfit``, ``factor1``, ``dist1`` ...).
    """
    fitSettings = resolve_fit_settings(fit_settings)
    rng = get_rng(rng)
    flag1 = flagX
    npatches_fixed = int(npatches_fixed)
    factor1 = None
    odeOptions = subepidemicOdeOptions()

    # <==============================================================================>
    # <============ Load data and proceed to parameter estimation ===================>
    # <==============================================================================>
    data1 = _prepare_window(data1, epidemic_period)
    data = data1[:, 1].copy()

    # <==============================================================================>
    # <============ Set time vector (timevect) and initial condition (I0) ===========>
    # <==============================================================================>
    timevect = data1[:, 0].copy()
    I0 = float(data[0])  # initial condition

    # <==============================================================================>
    # <================= Set range of C_thr values (onset_thrs) =====================>
    # <==============================================================================>
    # equal-width discretization of C_thr
    onset_thrs = onsetThresholdGrid(data, smoothfactor1)

    # <=============================================================================>
    # <===== Set range of the possible number of subepidemics (1:npatches_fixed)====>
    # <=============================================================================>
    npatchess = list(range(1, npatches_fixed + 1))
    if onset_fixed == 1 or (len(npatchess) == 1 and npatchess[0] == 1):
        onset_thrs = np.array([0.0])
    onset_thrs2 = onset_thrs

    nzMax = npatches_fixed * 4 + 2
    rankRows = []  # [nCandidate onset_thr AICc]
    paramRows = []

    # Diagnostics are in search order, including skipped partial-activation fits.
    thresholdDiagnostics = []
    fitDiagnostics = []
    candidateFitInfo = []
    candidatesTooLarge = []  # [sub-epidemics, parameters, data points] of sizes left out

    # <====================================================================================>
    # <==== Evaluate AICc across models with different number of subepidemics and C_thr ==>
    # <====================================================================================>
    # Series the objective is fitted to. Count likelihoods (methods 1, 3, 4, 5) are
    # fitted to the raw counts; least squares and absolute deviations to the moving
    # average. smoothfactor1 still shapes the threshold grid for every method.
    ydata, fitSeriesInfo = subepidemicFitSeries(data, smoothfactor1, method1, fitSettings)
    if fitSeriesInfo.smoothedCountLikelihood:
        # Only reached when countLikelihoodOnRawData = False in the fit settings.
        subepidemicWarning(
            "SubEpiPredict:SmoothedCountLikelihood",
            "method1 = %d fits a count likelihood to a %g-point moving average. The fitted dispersion "
            "will be too small and prediction intervals too narrow for the raw counts.",
            method1, smoothfactor1,
        )  # fmt: skip
    elif method1 in (1, 3, 4, 5) and smoothfactor1 > 1 and verbose:
        print(
            "method1 = %d is fitted to the raw counts. smoothfactor1 = %g is used only for the grid of "
            "onset thresholds." % (method1, smoothfactor1)
        )

    for nCandidate in npatchess:
        # The optimizer dimension is fixed for this entire threshold sweep.
        # A candidate that is too large for this series is not fitted: its AICc
        # would be infinite or, worse, lowered by a negative small-sample term.
        kCandidate = get_nparams(method1, dist1, nCandidate, flag1, 1, onset_fixed)
        if ydata.size <= kCandidate + 1:
            candidatesTooLarge.append([nCandidate, kCandidate, ydata.size])
            subepidemicWarning(
                "SubEpiPredict:CandidateTooLarge",
                "Models with %d sub-epidemic(s) have %d parameters, too many for the %d data points "
                "(AICc needs more than %d). They are left out of the screen.",
                nCandidate, kCandidate, ydata.size, kCandidate + 1,
            )  # fmt: skip
            continue

        if onset_fixed == 1 or nCandidate == 1:
            onset_thrs = np.array([0.0])
        else:
            onset_thrs = onset_thrs2.copy()
            # A threshold that the preceding sub-epidemic already meets at time
            # zero starts every sub-epidemic at once, exactly as threshold 0
            # does. Fitting it would rank the same model twice.
            levelAtStart = I0 if nCandidate == 2 else min(I0, 1)
            metAtStart = np.flatnonzero(onset_thrs <= levelAtStart)
            if metAtStart.size > 1:
                # All of these are one model; keep the first (threshold 0).
                sameAsZero = metAtStart[1:]
                if verbose:
                    print(
                        "Onset thresholds dropped as duplicates of threshold %g (already met at time zero): %s"
                        % (onset_thrs[metAtStart[0]], mat2str(onset_thrs[sameAsZero], 6))
                    )
                onset_thrs = np.delete(onset_thrs, sameAsZero)

        # <================================================================================================>
        # <=========================== Search space, seed and start points ===============================>
        # <================================================================================================>
        # Legacy bounds and seed unless the fit settings override them.
        LB, UB, isLog, z0 = subepidemicSearchSpace(nCandidate, data, flag1, method1, fitSettings)
        expectedDimension = 4 * nCandidate + 2

        # Start points ONCE per nCandidate, reused for every onset threshold:
        # numstartpoints Latin-hypercube points in unit coordinates (log-uniform
        # for r, K and alpha) and the seed. Rows are in original units.
        unitMap = subepidemicUnitMap("make", dict(z0=z0, LB=LB, UB=UB, isLog=isLog))
        nStarts = max(0, int(math.floor(numstartpoints)))
        starts_base = subepidemicStartPoints(unitMap, nStarts, 0, z0, rng)

        # Accepted solution of the previous threshold for this nCandidate
        # (continuation start); None until one candidate has been retained.
        previousP = None

        for onset_thr in onset_thrs:
            onset_thr = float(onset_thr)
            # Explicit fit context: the objective reads nothing from globals.
            context = Struct(
                timevect=timevect, ydata=ydata, I0=I0, npatches=nCandidate, onset_fixed=onset_fixed,
                onset_thr=onset_thr, flag1=flag1, method1=method1, odeOptions=odeOptions,
            )  # fmt: skip

            # The solution of the previous threshold is an excellent start for
            # this one. It is one more start point, never a replacement.
            starts = starts_base
            nStartsThisFit = nStarts
            if fitSettings.continuationStart and previousP is not None:
                starts = np.vstack([starts_base, previousP.reshape(1, -1)])
            elif previousP is None and fitSettings.firstThresholdStartFactor > 1 and nStarts > 0:
                # A fit without a solution of a previous threshold to start from
                # depends on the random start points alone. It gets
                # firstThresholdStartFactor times the start points (default 3).
                nExtra = int(math.ceil((fitSettings.firstThresholdStartFactor - 1) * nStarts))
                extraStarts = subepidemicStartPoints(unitMap, nExtra, 0, z0, rng)
                starts = np.vstack([extraStarts[:-1, :], starts_base])  # the seed is already the last row of starts_base
                nStartsThisFit = nStarts + nExtra

            P, fval, fitInfo = fitSubepidemicModel(
                context, z0, LB, UB, isLog, starts, nStartsThisFit, fitSettings, rng
            )

            candidateFitInfo.append(fitInfo)
            fitDiagnostics.append(
                [nCandidate, onset_thr, float(fitInfo.success), fitInfo.exitflag, fitInfo.startAgreement,
                 fitInfo.gradientInfNorm, float(fitInfo.polish.run), float(fitInfo.polish.improved),
                 fitInfo.nLocalRuns, fitInfo.nConvergedRuns,
                 np.size(fitInfo.atLowerBound) + np.size(fitInfo.atUpperBound)]
            )  # fmt: skip

            # A candidate with no converged, valid local minimum is recorded and
            # skipped. It is never replaced by its seed.
            if not fitInfo.success:
                thresholdDiagnostics.append([nCandidate, onset_thr, np.nan, 0, np.nan])
                subepidemicWarning(
                    "SubEpiPredict:CandidateRejected",
                    "Rejected n=%g, onset_thr=%g: %s", nCandidate, onset_thr, fitInfo.message,
                )  # fmt: skip
                continue

            # Activation is read from the accepted fit's own simulation.
            nActive = fitInfo.nActive
            thresholdDiagnostics.append([nCandidate, onset_thr, nActive, float(nActive == nCandidate), fval])

            # Conservative policy: do not relabel or rescore a partially active
            # nCandidate fit as a smaller model.
            if nActive != nCandidate:
                continue

            # Sub-epidemics that start together are interchangeable; give them a
            # fixed order. The curve and the AICc do not change.
            P, _, _ = subepidemicOrderLabels(P, nCandidate, onset_fixed, onset_thr, I0)
            previousP = P

            AICc = getAICc(method1, dist1, nCandidate, flag1, 1, fval, ydata.size, onset_fixed)[0]

            rankRows.append([nCandidate, onset_thr, AICc])
            row = np.zeros(nzMax)
            row[:expectedDimension] = P
            paramRows.append(row)

            if verbose:
                print(
                    "  %d sub-epidemic(s), C_thr = %-12s AICc = %.4f"
                    % (nCandidate, num2str(onset_thr), AICc)
                )

    thresholdDiagnostics = np.array(thresholdDiagnostics, dtype=float).reshape(-1, 5)
    fitDiagnostics = np.array(fitDiagnostics, dtype=float).reshape(-1, len(FIT_DIAGNOSTICS_COLUMNS))
    candidatesTooLarge = np.array(candidatesTooLarge, dtype=float).reshape(-1, 3)

    if not rankRows and candidatesTooLarge.shape[0] > 0 and fitDiagnostics.shape[0] == 0:
        raise SubEpiPredictError(
            "SubEpiPredict:TooFewDataPoints",
            "No model can be ranked: the series has %d data points and the smallest model has %d parameters. "
            "AICc needs more data points than parameters plus one."
            % (ydata.size, int(np.min(candidatesTooLarge[:, 1]))),
        )
    if not rankRows:
        raise SubEpiPredictError(
            "SubEpiPredict:NoFullyActiveCandidates",
            "No candidate activated all of its fitted subepidemics during calibration.",
        )

    # <=============================================================================================>
    # <======================== Sort the results by AICc (lowest to highest) =======================>
    # <=============================================================================================>
    RMSES = np.array(rankRows, dtype=float)
    PS = np.array(paramRows, dtype=float)
    index1 = np.lexsort((RMSES[:, 0], RMSES[:, 2]))  # sortrows(RMSES,[3 1])
    RMSES = RMSES[index1, :]
    PS = PS[index1, :]

    npatches = int(RMSES[0, 0])
    onset_thr = float(RMSES[0, 1])
    AICc_best = float(RMSES[0, 2])

    # The quantity exp((AICmin - AICi)/2) is the relative likelihood of model i.
    AICmin = RMSES[0, 2]
    relativelik_i = np.exp((AICmin - RMSES[:, 2]) / 2)

    # <=============================================================================================>
    # <============================ Get the best fit results =======================================>
    # <=============================================================================================>
    P = PS[0, : npatches * 4 + 2].copy()
    rs_hat, ps_hat, as_hat, Ks_hat, alpha_hat, d_hat = unpackSubepidemicParameters(P, npatches)

    if method1 in (3, 4, 5):
        dist1 = method1
        factor1 = alpha_hat

    _, x, totinc, patchIncidence, bestActivation = simulateSubepidemic(
        timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs_hat, ps_hat, as_hat, Ks_hat, odeOptions
    )

    # The stored dimension, parameter vector, and AICc must stay paired.
    nActiveBest = int(np.sum(bestActivation.active))
    if nActiveBest != npatches:
        raise SubEpiPredictError(
            "SubEpiPredict:BestModelActivationMismatch",
            "The selected %d-subepidemic model reconstructed with %d active patches. "
            "The saved model must be refitted, not truncated after ranking." % (npatches, nActiveBest),
        )

    bestfit = totinc

    # <=============================================================================================>
    # <================================= Plot best model fit ======================================>
    # <=============================================================================================>
    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure(100)
        fig.clf()
        ax = fig.add_subplot(1, 1, 1)
        for j in range(npatches):
            ax.plot(timevect, patchIncidence[:, j])
        ax.plot(timevect, totinc, "r")
        ax.plot(timevect, data, "bo", markerfacecolor="none")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Incidence")
        ax.set_title("best fit")

    if verbose:
        print("Best model of the screen: %d sub-epidemic(s), onset threshold %s." % (npatches, num2str(onset_thr)))

    if method1 == 0 and dist1 == 2:  # calculate the overdispersion factor
        # estimate dispersion in data
        binsize1 = 7
        ratios, _ = getMeanVarianceRatio(data, binsize1, 2)
        index1 = ratios[:, 0] > 0
        factor1 = float(np.mean(ratios[index1, 0])) if np.any(index1) else math.nan

    screen = Struct(
        thresholdDiagnostics=thresholdDiagnostics, thresholdDiagnosticsColumns=THRESHOLD_DIAGNOSTICS_COLUMNS,
        fitDiagnostics=fitDiagnostics, fitDiagnosticsColumns=FIT_DIAGNOSTICS_COLUMNS,
        candidateFitInfo=candidateFitInfo, candidatesTooLarge=candidatesTooLarge,
        relativelik_i=relativelik_i, AICc_best=AICc_best, timevect=timevect, I0=I0, data=data, ydata=ydata,
        data1=data1, onset_thrs2=onset_thrs2, bestfit=bestfit, factor1=factor1, dist1=dist1, method1=method1,
        flag1=flag1, npatches_fixed=npatches_fixed, onset_fixed=onset_fixed, smoothfactor1=smoothfactor1,
        calibrationperiod1=calibrationperiod1, fitSeriesInfo=fitSeriesInfo, file=None,
    )  # fmt: skip

    # <=============================================================================================>
    # <===================================  Save the results  ======================================>
    # <=============================================================================================>
    if save:
        from .io import save_mat

        os.makedirs(output_dir, exist_ok=True)
        name = "ABC-ensem-npatchesfixed-%s-onsetfixed-%s-smoothing-%s-%s-flag1-%s-method-%s-dist-%s-calibrationperiod-%s.mat" % (
            num2str(npatches_fixed), num2str(onset_fixed), num2str(smoothfactor1), str(datafilename1)[:-4],
            num2str(flag1), num2str(method1), num2str(dist1), num2str(calibrationperiod1),
        )  # fmt: skip
        screen.file = os.path.join(output_dir, name)
        save_mat(
            screen.file,
            dict(
                RMSES=RMSES, PS=PS, npatches=npatches, onset_thr=onset_thr, P=P, AICc_best=AICc_best,
                relativelik_i=relativelik_i.reshape(-1, 1), timevect=timevect.reshape(-1, 1), I0=I0,
                data=data.reshape(-1, 1), ydata=ydata.reshape(-1, 1), data1=data1,
                bestfit=bestfit.reshape(-1, 1), thresholdDiagnostics=thresholdDiagnostics,
                thresholdDiagnosticsColumns=THRESHOLD_DIAGNOSTICS_COLUMNS, fitDiagnostics=fitDiagnostics,
                fitDiagnosticsColumns=FIT_DIAGNOSTICS_COLUMNS, candidatesTooLarge=candidatesTooLarge,
                onset_thrs2=onset_thrs2, factor1=factor1, dist1=dist1, method1=method1, flag1=flag1,
                npatches_fixed=npatches_fixed, onset_fixed=onset_fixed, smoothfactor1=smoothfactor1,
                calibrationperiod1=calibrationperiod1, numstartpoints=numstartpoints, M=M, DT=DT,
                datafilename1=str(datafilename1),
            ),
        )  # fmt: skip

    return RMSES, PS, npatches, onset_thr, P, screen


# <============================================================================>
# <=============================== Bootstrap ==================================>
# <============================================================================>


def fittingModifiedLogisticFunctionPatchMultiple(
    RMSES, PS, data1, DT, epidemic_period, M, flagX, numstartpoints, rank1, *,
    method1, dist1, onset_fixed, smoothfactor1, factor1=None,
    fit_settings=None, plot=True, verbose=True, rng=None,
):  # fmt: skip
    """Parametric bootstrap of one ranked model.

    Port of ``fittingModifiedLogisticFunctionPatchMultiple.m``. The ranked fit
    is bootstrapped as it stands: ``M`` synthetic datasets are drawn around
    its curve with the error structure ``dist1`` (:func:`AddPoissonError`) and
    each is refitted by :func:`fitSubepidemicModel` from the ranked vector,
    one jittered copy of it and ``bootstrapStarts`` Latin-hypercube points. A
    refit with no converged local minimum is counted as unresolved and its
    best valid end point is kept.

    Parameters
    ----------
    RMSES, PS : ranking and parameter vectors of the screen
    data1 : n-by-2 array, time index and incidence
    DT, numstartpoints : kept as in MATLAB (not used)
    epidemic_period : ONE-BASED rows of ``data1`` to use
    M : number of bootstrap realizations
    flagX : growth model
    rank1 : rank of the model to bootstrap (ONE-BASED: 1 is the best model)
    method1, dist1, onset_fixed, smoothfactor1 : the settings MATLAB reads
        from global variables
    factor1 : variance-to-mean factor for ``method1 = 0`` with ``dist1 = 2``
        (estimated by the driver); for the other error structures it is
        computed here

    Returns
    -------
    ``(Phatss, npatches, onset_thr, curves, bestfit, data1, P0, AICc_best, factor1, d, bootstrapDiagnostics)``
        ``Phatss``: M-by-(4n+2) refitted vectors. ``curves``: the synthetic
        datasets (columns). ``bestfit``: the curve of the ranked fit. ``P0``:
        the ranked vector. ``factor1``: the noise scale used. ``d``: the power
        d of the ranked fit. ``bootstrapDiagnostics``: Struct with one entry
        per refit (``success``, ``usedFallback``, ``exitflag``, ``nActive``,
        ``startAgreement``, ``gradientInfNorm``, ``polishRun``,
        ``polishImproved``, ``nAtBound``, ``fval``) and ``rankedAtBound``,
        ``boundReport``, ``searchLB``, ``searchUB``, ``searchIsLog``,
        ``identifiability``, ``notIdentified``.
    """
    fitSettings = resolve_fit_settings(fit_settings)
    rng = get_rng(rng)
    flag1 = flagX
    M = int(M)
    rank1 = int(rank1)

    # remove repeated rows
    RMSES = np.asarray(RMSES, dtype=float)
    PS = np.asarray(PS, dtype=float)
    RMSES, index1 = unique_rows_stable(RMSES)
    PS = PS[index1, :]

    npatches = int(RMSES[rank1 - 1, 0])
    onset_thr = float(RMSES[rank1 - 1, 1])
    AICc_best = float(RMSES[rank1 - 1, 2])
    P0 = PS[rank1 - 1, : npatches * 4 + 2].copy()

    numparams = get_nparams(method1, dist1, npatches, flag1, 1, onset_fixed)

    data1 = _prepare_window(data1, epidemic_period)
    data = data1[:, 1].copy()
    timevect = data1[:, 0].copy()

    # The ranked vector P0 is used unchanged and its score is checked.
    d = float(P0[-1])  # legacy output: the power d of the ranked fit
    odeOptions = subepidemicOdeOptions()

    # Same bounds and coordinate scales as the screen used for this candidate.
    LB, UB, isLog, _ = subepidemicSearchSpace(npatches, data, flag1, method1, fitSettings)

    I0 = float(data[0])  # initial condition

    # Same series as the screen fitted (raw counts for count likelihoods).
    ydata, _ = subepidemicFitSeries(data, smoothfactor1, method1, fitSettings)

    # Explicit fit context: the objective reads nothing from globals.
    context = Struct(
        timevect=timevect, ydata=ydata, I0=I0, npatches=npatches, onset_fixed=onset_fixed, onset_thr=onset_thr,
        flag1=flag1, method1=method1, odeOptions=odeOptions,
    )  # fmt: skip
    fvalRanked, rankedDetail = subepidemicObjective(P0, context)
    if not rankedDetail.success:
        raise SubEpiPredictError(
            "SubEpiPredict:RankedFitInvalid",
            "The ranked model %d could not be reproduced: %s" % (rank1, rankedDetail.message),
        )
    if rankedDetail.nActive != npatches:
        raise SubEpiPredictError(
            "SubEpiPredict:BestModelActivationMismatch",
            "Ranked model %d was fitted with %d sub-epidemics but reconstructs with %d active."
            % (rank1, npatches, rankedDetail.nActive),
        )
    AICcCheck = getAICc(method1, dist1, npatches, flag1, 1, fvalRanked, ydata.size, onset_fixed)[0]
    if abs(AICcCheck - AICc_best) > 1e-8 * max(1, abs(AICc_best)):
        subepidemicWarning(
            "SubEpiPredict:RankedScoreChanged",
            "Ranked model %d: AICc %.10g at the saved vector differs from the ranked AICc %.10g. "
            "Data, smoothing or ODE settings changed between the screen and the bootstrap.",
            rank1, AICcCheck, AICc_best,
        )  # fmt: skip

    # A ranked fit with r or K at a search limit is set by the limit, not by the data.
    rankedAtBound, _, _ = subepidemicBoundSummary(P0, LB, UB, isLog, npatches)
    for item in rankedAtBound:
        if item.isSearchLimit:
            subepidemicWarning(
                "SubEpiPredict:RankedModelAtBound",
                "Ranked model %d has %s at its %s limit (%g). The estimate is set by the limit. "
                "Set %s in the fit settings to move it.",
                rank1, item.name, item.side, item.limit, item.setting,
            )  # fmt: skip

    P = P0
    rs_hat = P[0:npatches]
    ps_hat = P[npatches : 2 * npatches]
    as_hat = P[2 * npatches : 3 * npatches]
    Ks_hat = P[3 * npatches : 4 * npatches]
    alpha_hat = float(P[-2])
    d_hat = float(P[-1])

    _, x, totinc, patchIncidence, _ = simulateSubepidemic(
        timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs_hat, ps_hat, as_hat, Ks_hat, odeOptions
    )

    bestfit = totinc

    if method1 == 0 and dist1 == 0:  # normal distribution of the error structure
        var1 = np.sum((bestfit - data) ** 2) / (bestfit.size - numparams)
        factor1 = float(np.sqrt(var1))
    elif method1 == 6 and dist1 == 6:
        var1 = np.sum((bestfit - data) ** 2) / (bestfit.size - numparams)
        factor1 = float(np.sqrt(var1) / np.sqrt(2))  # scale parameter: b = sigma/sqrt(2)
    elif method1 in (3, 4, 5):
        # The synthetic datasets of this ranked model are drawn with ITS OWN
        # fitted dispersion.
        factor1 = alpha_hat

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure(10)
        fig.clf()
        ax = fig.add_subplot(1, 1, 1)
        for j in range(npatches):
            ax.plot(timevect, patchIncidence[:, j])
        ax.plot(timevect, totinc, "r", label="Sub-epidemics=%s; $C_{thr}$=%s" % (num2str(npatches), num2str(onset_thr)))
        ax.plot(timevect, data, "bo", markerfacecolor="none")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("incidence fit")
        ax.set_title("Num. Subepidemics=%s; AICc=%s" % (num2str(npatches), num2str(AICc_best, 6)))
        ax.legend()

    # <=============================================================================================>
    # <======================== Vector with parameter estimates ====================================>
    # <=============================================================================================>
    Ptrue = np.concatenate([rs_hat, ps_hat, as_hat, Ks_hat, [alpha_hat, d_hat]])
    if verbose:
        print("\nRanked model %d, estimates [r p a K alpha d]: %s" % (rank1, mat2str(Ptrue, 5)))
        print("Bootstrap: %d refits of ranked model %d ..." % (M, rank1))

    # <=============================================================================================>
    # <========== Generate parameter uncertainty via parametric bootstrapping ======================>
    # <=============================================================================================>
    yi = np.cumsum(totinc)
    Phatss = np.zeros((M, 4 * npatches + 2))
    curves = np.zeros((yi.size, M))

    bootstrapDiagnostics = Struct(
        success=np.zeros(M, dtype=bool), usedFallback=np.zeros(M, dtype=bool), exitflag=np.full(M, np.nan),
        nActive=np.full(M, np.nan), startAgreement=np.full(M, np.nan), gradientInfNorm=np.full(M, np.nan),
        polishRun=np.zeros(M, dtype=bool), polishImproved=np.zeros(M, dtype=bool), nAtBound=np.full(M, np.nan),
        fval=np.full(M, np.nan),
    )  # fmt: skip
    bootContext = Struct(context)

    # Bootstrap counts for user-facing diagnostics only.
    acceptedBootstrap = 0

    for real in range(M):
        yirData = AddPoissonError(yi, 1, dist1, factor1, d_hat, rng)[:, 0]
        curves[:, real] = yirData
        bootContext.ydata = yirData

        Pboot, fvalBoot, bootInfo = fitSubepidemicModel(
            bootContext, P, LB, UB, isLog, None, fitSettings.bootstrapStarts, fitSettings, rng
        )
        if bootInfo.success:
            nActiveBoot = bootInfo.nActive
        elif bootInfo.fallback is not None:
            # unresolved refit: keep its best valid end point (legacy policy)
            Pboot = bootInfo.fallback.P
            fvalBoot = bootInfo.fallback.fval
            nActiveBoot = bootInfo.fallback.nActive
            bootstrapDiagnostics.usedFallback[real] = True
        else:
            raise SubEpiPredictError(
                "SubEpiPredict:BootstrapFitFailed",
                "Bootstrap dataset %d of ranked model %d has no valid parameter vector: %s"
                % (real + 1, rank1, bootInfo.message),
            )

        acceptedBootstrap += int(bootInfo.success)
        bootstrapDiagnostics.success[real] = bootInfo.success
        bootstrapDiagnostics.exitflag[real] = bootInfo.exitflag
        bootstrapDiagnostics.nActive[real] = nActiveBoot
        bootstrapDiagnostics.startAgreement[real] = bootInfo.startAgreement
        bootstrapDiagnostics.gradientInfNorm[real] = bootInfo.gradientInfNorm
        bootstrapDiagnostics.polishRun[real] = bootInfo.polish.run
        bootstrapDiagnostics.polishImproved[real] = bootInfo.polish.improved
        bootstrapDiagnostics.nAtBound[real] = np.size(bootInfo.atLowerBound) + np.size(bootInfo.atUpperBound)
        bootstrapDiagnostics.fval[real] = fvalBoot

        # The same fixed order of interchangeable sub-epidemics as in the
        # screen, so a refit cannot swap their labels.
        Phatss[real, :], _, _ = subepidemicOrderLabels(Pboot, npatches, onset_fixed, onset_thr, data[0])

    boundReport, boundLines, nRefitsAtBound = subepidemicBoundSummary(Phatss, LB, UB, isLog, npatches)
    if verbose:
        print("\nBootstrap summary for ranked model %d:" % rank1)
        print("Requested bootstrap datasets: %d" % M)
        print("Accepted bootstrap datasets: %d" % acceptedBootstrap)
        print("Unresolved bootstrap datasets: %d" % (M - acceptedBootstrap))
        print(
            "Refits with fewer than %d active sub-epidemics: %d"
            % (npatches, int(np.sum(bootstrapDiagnostics.nActive < npatches)))
        )
        print("Refits with a parameter at a bound: %d" % nRefitsAtBound)
        for line in boundLines:
            print("   %s refits" % line)
        print("Refits in which the polish ran: %d\n" % int(np.sum(bootstrapDiagnostics.polishRun)))
    bootstrapDiagnostics.rankedAtBound = rankedAtBound
    bootstrapDiagnostics.boundReport = boundReport

    # <================================================================================================>
    # <=========================== Parameter identifiability (PII) ====================================>
    # <================================================================================================>
    bootstrapDiagnostics.searchLB = LB
    bootstrapDiagnostics.searchUB = UB
    bootstrapDiagnostics.searchIsLog = isLog
    bootstrapDiagnostics.identifiability = subepidemicIdentifiability(
        Phatss, LB, UB, npatches, P0, fitSettings.PIIscale
    )
    if verbose:
        subepidemicPrintPII(
            bootstrapDiagnostics.identifiability, "Parameter identifiability (PII) of ranked model %d" % rank1
        )

    # A parameter whose 95% interval covers at least notIdentifiedRatio of what
    # the search bounds allow is not identified by the data.
    piiReport = bootstrapDiagnostics.identifiability
    with np.errstate(invalid="ignore"):
        notIdentified = np.flatnonzero(~piiReport.fixed & (piiReport.PII_ratio >= fitSettings.notIdentifiedRatio))
    for k in notIdentified:
        if piiReport.PII_scale[k] == 1:
            coverText = "%.0f%% of the search range on the log scale (%.1f of %.1f decades)" % (
                100 * piiReport.PII_ratio[k], piiReport.PII[k], piiReport.PII_box[k],
            )  # fmt: skip
        else:
            coverText = "%.0f%% of the search range" % (100 * piiReport.PII_ratio[k])
        subepidemicWarning(
            "SubEpiPredict:NotIdentified",
            "Ranked model %d: %s is not identified by these data. Its 95%% interval, %.4g to %.4g, covers %s, "
            "so it reflects the search bounds, not the data; moving a bound moves the interval. "
            "Do not interpret its estimate (%.4g).",
            rank1, piiReport.name[k], piiReport.CI_LB[k], piiReport.CI_UB[k], coverText, piiReport.estimate[k],
        )  # fmt: skip
    bootstrapDiagnostics.notIdentified = [piiReport.name[k] for k in notIdentified]
    notIdentifiedColumns = [int(piiReport.column[k]) for k in notIdentified]

    # More than 2.5% of the refits at a limit of r or K: the 95% percentile
    # interval of that parameter ends at the limit instead of where the data put it.
    for item in boundReport:
        if item.isSearchLimit and item.fraction > 0.025 and item.index not in notIdentifiedColumns:
            subepidemicWarning(
                "SubEpiPredict:BootstrapAtBound",
                "Ranked model %d: %d of %d bootstrap refits have %s at its %s limit (%g), so its 95%% interval "
                "ends at the limit. Set %s in the fit settings to move it.",
                rank1, item.count, item.total, item.name, item.side, item.limit, item.setting,
            )  # fmt: skip

    return Phatss, npatches, onset_thr, curves, bestfit, data1, P0, AICc_best, factor1, d, bootstrapDiagnostics
