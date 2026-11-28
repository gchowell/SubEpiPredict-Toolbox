"""The n-sub-epidemic model: growth kernels, initial state, event-driven simulation.

Ports of ``modifiedLogisticGrowthPatch.m``, ``initializeSubepidemicState.m``,
``unpackSubepidemicParameters.m``, ``subepidemicIncidence.m``,
``subepidemicOdeOptions.m``, ``subepidemicActivationEvents.m``,
``simulateSubepidemic.m`` and ``subepidemicSimulateRefit.m``.

How the simulation is computed here
-----------------------------------
MATLAB integrates all sub-epidemics together with ``ode15s`` and stops at the
first upward crossing of the onset threshold (an ODE event), switches the next
sub-epidemic on and restarts. The equations of the sub-epidemics are not
coupled: sub-epidemic ``j`` only needs to know WHEN sub-epidemic ``j-1``
reached the threshold. The port therefore integrates one sub-epidemic at a
time with LSODA (``scipy.integrate.odeint``), finds the crossing time of the
threshold to the solver tolerance, and starts the next sub-epidemic from its
seed at that time. The trajectory is the same function of time as in MATLAB;
the two solvers agree to their tolerances (RelTol 1e-8, AbsTol 1e-10).

One expression is guarded. The generalized Richards kernel (flag1 = 2) contains
``(1 - C/K)^a``, which is complex when ``C > K`` and ``a`` is not a whole
number. The solution of the model never exceeds ``K``, but a solver probes
states slightly above it once a sub-epidemic has saturated. The port sets the
growth rate to zero there, so a saturated curve stays at ``K``; MATLAB returns
a complex trajectory in that case, which ``simulateSubepidemic.m`` refuses.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

from .errors import SubEpiPredictError
from .utils import Struct

_ODE_OK = "Integration successful."

# Identifiers of the failures that mean "this parameter set cannot be
# simulated" (subepidemicSimulateRefit.m).
INTEGRATION_FAILURES = (
    "SubEpiPredict:IncompleteIntegration",
    "SubEpiPredict:InvalidTrajectory",
    "SubEpiPredict:InvalidSegment",
    "SubEpiPredict:MissedActivationEvent",
    "SubEpiPredict:ActivationProgressFailure",
)


# <============================================================================>
# <============================== Small validators ============================>
# <============================================================================>


def _is_real_number(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def _check_patch_count(npatches, where=""):
    if not (_is_real_number(npatches) and math.isfinite(npatches) and npatches >= 1 and npatches == int(npatches)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidPatchCount", "npatches must be a positive integer%s." % where
        )
    return int(npatches)


def _check_onset_mode(onset_fixed):
    ok = isinstance(onset_fixed, (bool, np.bool_)) or _is_real_number(onset_fixed)
    if not (ok and (onset_fixed == 0 or onset_fixed == 1)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidOnsetMode", "onset_fixed must be 0 (asynchronous) or 1 (synchronous)."
        )
    return int(onset_fixed)


def _check_flag(flag1):
    if not (_is_real_number(flag1) and math.isfinite(flag1) and flag1 == int(flag1) and 0 <= flag1 <= 5):
        raise SubEpiPredictError("SubEpiPredict:InvalidModelFlag", "flag1 must be an integer from 0 to 5.")
    return int(flag1)


# <============================================================================>
# <=============================== Growth kernels =============================>
# <============================================================================>


def _rhs(flag1, r, p, a, K):
    """Right-hand side ``rhs(y, t)`` of ONE sub-epidemic for ``odeint`` (y has one entry).

    The six expressions of ``modifiedLogisticGrowthPatch.m``. A value that
    MATLAB would return as NaN or complex (negative base with a fractional
    exponent, log of a negative number, division by zero) is returned as NaN,
    and an overflow as Inf; both make the simulation fail, as in MATLAB.
    """
    pow_ = math.pow
    log = math.log
    nan = math.nan
    inf = math.inf

    if flag1 == 0:  # GGM

        def rhs(y, t):
            try:
                return (r * pow_(y[0], p),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    elif flag1 == 1:  # GLM

        def rhs(y, t):
            c = y[0]
            try:
                return ((r * pow_(c, p)) * (1 - c / K),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    elif flag1 == 2:  # Generalized Richards (existing expression retained)
        aWhole = float(a).is_integer()

        def rhs(y, t):
            c = y[0]
            try:
                base = 1 - c / K
                if base < 0 and not aWhole:
                    # (negative)^a is complex for a fractional a. The solution of
                    # the model never exceeds K, but a solver may probe a state
                    # just above it; growth is zero there (see the module notes).
                    base = 0.0
                return ((r * pow_(c, p)) * pow_(base, a),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    elif flag1 == 3:  # Logistic

        def rhs(y, t):
            c = y[0]
            try:
                return ((r * c) * (1 - c / K),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    elif flag1 == 4:  # Richards

        def rhs(y, t):
            c = y[0]
            try:
                return ((r * c) * (1 - pow_(c / K, a)),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    elif flag1 == 5:  # Gompertz

        def rhs(y, t):
            c = y[0]
            try:
                return (r * c * log(K / c),)
            except OverflowError:
                return (inf,)
            except (ValueError, ZeroDivisionError):
                return (nan,)

    else:
        raise SubEpiPredictError("SubEpiPredict:InvalidModelFlag", "flag1 must be an integer from 0 to 5.")

    return rhs


def _kernel(flag1, r, p, a, K):
    """Growth rate dC/dt = f(C) of ONE sub-epidemic as a function of a float."""
    rhs = _rhs(flag1, r, p, a, K)
    return lambda c: rhs((c,), 0.0)[0]


def modifiedLogisticGrowthPatch(t, x, rs1, ps1, as1, Ks1, npatches, onset_thr, flag1, activeMask=None):
    """Derivative of every sub-epidemic for one fixed activation mask.

    Port of ``modifiedLogisticGrowthPatch.m``. ``activeMask`` holds one zero
    or one per sub-epidemic; a dormant sub-epidemic has derivative 0 and its
    growth expression is not evaluated. ``t`` and ``onset_thr`` remain in the
    interface but do not determine the activation: threshold decisions belong
    to :func:`simulateSubepidemic`. Returns an array of length ``npatches``.

    Growth expressions (flag1): 0 ``r C^p``; 1 ``r C^p (1 - C/K)``;
    2 ``r C^p (1 - C/K)^a``; 3 ``r C (1 - C/K)``; 4 ``r C (1 - (C/K)^a)``;
    5 ``r C log(K/C)``.
    """
    if activeMask is None:
        raise SubEpiPredictError(
            "SubEpiPredict:MissingActivationMask", "Supply an explicit activeMask, or call simulateSubepidemic."
        )
    mask = np.asarray(activeMask)
    ok = mask.ndim <= 1 and mask.size == npatches and mask.dtype.kind in "biuf"
    if ok:
        mask = mask.astype(float).ravel()
        ok = bool(np.all((mask == 0) | (mask == 1)))
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidActivationMask", "activeMask must contain one zero or one for each patch."
        )
    if flag1 not in (0, 1, 2, 3, 4, 5):
        raise SubEpiPredictError("SubEpiPredict:InvalidModelFlag", "flag1 must be an integer from 0 to 5.")
    x = np.asarray(x, dtype=float).ravel()
    rs1 = np.asarray(rs1, dtype=float).ravel()
    ps1 = np.asarray(ps1, dtype=float).ravel()
    as1 = np.asarray(as1, dtype=float).ravel()
    Ks1 = np.asarray(Ks1, dtype=float).ravel()
    dx = np.zeros(int(npatches))
    for j in range(int(npatches)):
        if not mask[j]:
            # Do not evaluate a dormant growth expression: 0*NaN is still NaN.
            continue
        dx[j] = _kernel(int(flag1), float(rs1[j]), float(ps1[j]), float(as1[j]), float(Ks1[j]))(float(x[j]))
    return dx


# <============================================================================>
# <=========================== State and parameters ===========================>
# <============================================================================>


def initializeSubepidemicState(I0, npatches, onset_fixed):
    """Initial state for every sub-epidemic solve.

    Port of ``initializeSubepidemicState.m``. ``onset_fixed``, not the
    threshold value, defines the initialization:

    0  asynchronous: ``IC = [I0, 1, ..., 1]``, initially only the first
       sub-epidemic is active;
    1  synchronous: ``IC = I0/npatches`` for all, all active.

    The one-case seeds of later asynchronous sub-epidemics are numerical
    initial states, not additional cases in the first incidence observation.

    Returns ``(IC, invasions0, timeinvasions0, Cinvasions0)``.
    """
    if not (_is_real_number(I0) and math.isfinite(I0) and I0 >= 0):
        raise SubEpiPredictError("SubEpiPredict:InvalidI0", "I0 must be a finite, real, nonnegative scalar.")
    npatches = _check_patch_count(npatches)
    onset_fixed = _check_onset_mode(onset_fixed)

    IC = np.ones(npatches)
    invasions0 = np.zeros(npatches)
    timeinvasions0 = np.zeros(npatches)
    Cinvasions0 = np.zeros(npatches)

    if onset_fixed == 0:
        IC[0] = I0
        invasions0[0] = 1
    else:
        IC[:] = I0 / npatches
        invasions0[:] = 1
    return IC, invasions0, timeinvasions0, Cinvasions0


def unpackSubepidemicParameters(z, npatches):
    """Decode a parameter vector ``[r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]``.

    Port of ``unpackSubepidemicParameters.m``. Returns
    ``(rs, ps, as_, Ks, alpha, d)``; the four growth blocks are arrays of
    length ``npatches``.
    """
    npatches = _check_patch_count(npatches, " when decoding parameters")
    expectedDimension = 4 * npatches + 2
    try:
        zz = np.asarray(z, dtype=float)
    except (TypeError, ValueError):
        zz = None
    if zz is None or zz.ndim > 2 or zz.size != expectedDimension or (zz.ndim == 2 and 1 not in zz.shape):
        received = 0 if zz is None else zz.size
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Expected a real parameter vector of length %d for %d subepidemics; received %d entries."
            % (expectedDimension, npatches, received),
        )
    zz = zz.ravel()
    n = npatches
    return zz[0:n], zz[n : 2 * n], zz[2 * n : 3 * n], zz[3 * n : 4 * n], float(zz[-2]), float(zz[-1])


def subepidemicIncidence(x, onset_fixed):
    """Convert a full raw ODE trajectory to interval incidence.

    Port of ``subepidemicIncidence.m``. ``x`` has one row per observation time
    and one column per sub-epidemic and must begin at the calibration initial
    time. Later rows are interval increments. Only the first row is corrected
    for the one-case seeds of the asynchronous sub-epidemics 2..n.

    Returns ``(totinc, patchIncidence)``.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    if x.ndim != 2 or x.size == 0:
        raise SubEpiPredictError("SubEpiPredict:InvalidTrajectory", "x must be a nonempty numeric matrix.")
    onset_fixed = _check_onset_mode(onset_fixed)
    patchIncidence = np.vstack([x[0:1, :], np.diff(x, axis=0)])
    if onset_fixed == 0:
        patchIncidence[0, 1:] = patchIncidence[0, 1:] - 1
    totinc = patchIncidence.sum(axis=1)
    return totinc, patchIncidence


def subepidemicOdeOptions():
    """One place for the ODE solver settings of the toolbox.

    Port of ``subepidemicOdeOptions.m``: ``RelTol = 1e-8``, ``AbsTol = 1e-10``
    for every fit, bootstrap refit, forecast and plot. Incidence is the
    difference of two cumulative values that are much larger, so the solver
    defaults (1e-3, 1e-6) give AICc values that are off by as much as the gaps
    between the top-ranked models, and noisy finite-difference gradients.
    """
    return Struct(RelTol=1e-8, AbsTol=1e-10)


def subepidemicActivationEvents(t, x, onset_thr, watchedPatches):
    """Upward threshold-crossing event function (``subepidemicActivationEvents.m``).

    ``watchedPatches`` holds the dormant sub-epidemics (numbered from 1, all
    >= 2). Returns ``(value, isterminal, direction)``; ``value[k]`` is the
    state of the predecessor of ``watchedPatches[k]`` minus the threshold.
    The simulation of the port does not call it (see the module notes); it is
    kept so that the definition of an activation is in one place.
    """
    x = np.asarray(x, dtype=float).ravel()
    watched = np.asarray(watchedPatches, dtype=int).ravel()
    value = x[watched - 2] - onset_thr
    return value, np.ones(watched.size), np.ones(watched.size)  # upward crossings only


# <============================================================================>
# <================================ Simulation ================================>
# <============================================================================>


def _integrate(rhs, c0, tgrid, rtol, atol):
    """Integrate dC/dt = rhs from ``tgrid[0]`` and return C at ``tgrid`` (array)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with np.errstate(all="ignore"):
            y, info = odeint(rhs, (c0,), tgrid, rtol=rtol, atol=atol, mxstep=5000, full_output=True)
    if info["message"] != _ODE_OK:
        raise SubEpiPredictError(
            "SubEpiPredict:IncompleteIntegration",
            "The ODE solver stopped before the requested endpoint without an activation event.",
        )
    return y[:, 0]


def _crossing_time(rhs, tA, cA, tB, thr, rtol, atol):
    """Time in (tA, tB] at which the solution through (tA, cA) reaches ``thr``.

    The solution of an autonomous scalar equation is monotone, so time is a
    function of the state there: dt/dC = 1/f(C), integrated from cA to thr.
    If that fails, the time is found by bisection on the forward solution.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with np.errstate(all="ignore"):
                y, info = odeint(
                    lambda t, c: (1.0 / rhs((c,), 0.0)[0],),
                    (tA,),
                    (cA, thr),
                    rtol=min(rtol, 1e-10),
                    atol=1e-12,
                    mxstep=5000,
                    full_output=True,
                )
        tau = float(y[1, 0])
        slack = 1e-9 * max(1.0, abs(tB))
        if info["message"] == _ODE_OK and math.isfinite(tau) and tA - slack <= tau <= tB + slack:
            return min(max(tau, tA), tB)
    except (ZeroDivisionError, FloatingPointError, ValueError):
        pass

    def g(t):
        if t <= tA:
            return cA - thr
        return _integrate(rhs, cA, (tA, t), rtol, atol)[1] - thr

    if g(tB) < 0:
        return tB
    return float(brentq(g, tA, tB, xtol=1e-12 * max(1.0, abs(tB)), rtol=8.9e-16, maxiter=200))


def _simulate_core(t, I0, n, onset_fixed, onset_thr, flag1, rs, ps, as_, Ks, rtol, atol):
    """Simulation without input checks. Returns ``(x, active, time, count, initialActive, nSegments)``."""
    nt = t.size
    t0 = t[0]
    x = np.empty((nt, n))
    active = np.zeros(n, dtype=bool)
    actTime = np.full(n, np.nan)
    actCount = np.full(n, np.nan)
    initialActive = np.zeros(n, dtype=bool)

    if onset_fixed == 1:
        c0 = I0 / n
        for j in range(n):
            f = _rhs(flag1, rs[j], ps[j], as_[j], Ks[j])
            x[:, j] = _integrate(f, c0, t, rtol, atol) if nt > 1 else c0
        active[:] = True
        initialActive[:] = True
        actTime[:] = t0
        nSegments = 1 if nt > 1 else 0
    else:
        initialActive[0] = True
        eventTimes = set()
        fPrev = None
        tauPrev = t0  # activation time of the predecessor (None: never activated)
        for j in range(n):
            f = _rhs(flag1, rs[j], ps[j], as_[j], Ks[j])
            if j == 0:
                tau = t0
                c0 = I0
            else:
                c0 = 1.0
                prev = x[:, j - 1]
                # Literal rule of the toolbox: the predecessor's STATE is compared
                # with the threshold, also while the predecessor is a dormant seed.
                if prev[0] >= onset_thr:
                    tau = t0
                    actCount[j] = prev[0]
                else:
                    hit = np.flatnonzero(prev >= onset_thr)
                    if hit.size == 0 or tauPrev is None:
                        tau = None
                    else:
                        k = int(hit[0])
                        tA = t[k - 1]
                        cA = prev[k - 1]
                        if tauPrev > tA:  # the predecessor itself started inside this interval
                            tA = tauPrev
                            cA = 1.0 if j - 1 > 0 else I0
                        tau = _crossing_time(fPrev, tA, cA, t[k], onset_thr, rtol, atol)
                        actCount[j] = onset_thr
                        if tau < t[-1]:
                            eventTimes.add(tau)
            if tau is None:
                x[:, j] = c0
            else:
                active[j] = True
                actTime[j] = tau
                later = t > tau
                x[~later, j] = c0
                if np.any(later):
                    grid = np.concatenate(([tau], t[later]))
                    x[later, j] = _integrate(f, c0, grid, rtol, atol)[1:]
            fPrev = f
            tauPrev = tau
        nSegments = (1 + len(eventTimes)) if nt > 1 else 0

    if not np.all(np.isfinite(x)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidTrajectory",
            "Simulation must return a finite, real state at every requested time.",
        )
    return x, active, actTime, actCount, initialActive, nSegments


def simulateSubepidemic(
    timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs1, ps1, as1, Ks1, odeOptions=None
):
    """Simulate the n-sub-epidemic model on a time grid.

    Port of ``simulateSubepidemic.m``, shared by the fit and the forecast.

    Parameters
    ----------
    timevect : strictly increasing times; the first one is the calibration
        initial time (not a forecast-only tail)
    I0 : initial value (first observation)
    npatches : number of sub-epidemics n
    onset_fixed : 0 asynchronous (sub-epidemic j starts when sub-epidemic j-1
        reaches ``onset_thr``), 1 synchronous (all start at once from ``I0/n``)
    onset_thr : onset threshold C_thr (cumulative count of the predecessor)
    flag1 : growth model, 0 to 5 (see :func:`modifiedLogisticGrowthPatch`)
    rs1, ps1, as1, Ks1 : one value per sub-epidemic
    odeOptions : :func:`subepidemicOdeOptions` when omitted; only ``RelTol``
        and ``AbsTol`` are used. ``Events``, ``Mass``, ``Jacobian`` and
        ``Vectorized`` are refused, as in MATLAB.

    Returns
    -------
    ``(t, x, totinc, patchIncidence, activation)``
        ``t`` the grid; ``x`` the cumulative curve of every sub-epidemic
        (rows = times); ``totinc`` the total incidence per interval;
        ``patchIncidence`` the incidence of every sub-epidemic; ``activation``
        a Struct with ``active`` (final mask), ``time`` (initial time for
        initially active sub-epidemics, crossing time for triggered ones, NaN
        for never activated), ``predecessorCount`` (predecessor state at
        triggering), ``initialActive`` and ``nSegments``.

    A zero threshold does not turn an asynchronous model into a synchronous
    one: the first sub-epidemic still starts at ``I0`` and the others at 1.
    Event times are not inserted into the grid.
    """
    if odeOptions is None or (isinstance(odeOptions, dict) and len(odeOptions) == 0):
        odeOptions = subepidemicOdeOptions()

    n = _check_patch_count(npatches)
    onset_fixed = _check_onset_mode(onset_fixed)
    if not (_is_real_number(I0) and math.isfinite(I0) and I0 >= 0):
        raise SubEpiPredictError("SubEpiPredict:InvalidI0", "I0 must be a finite, real, nonnegative scalar.")

    try:
        t = np.asarray(timevect, dtype=float)
    except (TypeError, ValueError):
        t = np.empty(0)
    if t.ndim == 2 and 1 in t.shape:
        t = t.ravel()
    if t.ndim == 0:
        t = t.reshape(1)
    if not (t.ndim == 1 and t.size > 0 and np.all(np.isfinite(t)) and np.all(np.diff(t) > 0)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidTimeGrid", "timevect must be a finite, real, strictly increasing time vector."
        )
    if not (_is_real_number(onset_thr) and math.isfinite(onset_thr) and onset_thr >= 0):
        raise SubEpiPredictError("SubEpiPredict:InvalidThreshold", "onset_thr must be finite and nonnegative.")
    flag1 = _check_flag(flag1)

    parameters = []
    for v in (rs1, ps1, as1, Ks1):
        try:
            a = np.asarray(v, dtype=float)
        except (TypeError, ValueError):
            raise SubEpiPredictError(
                "SubEpiPredict:InvalidGrowthParameters", "Growth parameters must be finite, real numeric vectors."
            ) from None
        parameters.append(a)
    if any(a.size != n for a in parameters):
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Each growth-parameter vector must contain exactly npatches entries.",
        )
    for a in parameters:
        if not (a.ndim <= 2 and np.all(np.isfinite(a))):
            raise SubEpiPredictError(
                "SubEpiPredict:InvalidGrowthParameters", "Growth parameters must be finite, real numeric vectors."
            )
    rs, ps, as_, Ks = [[float(v) for v in a.ravel()] for a in parameters]

    if not isinstance(odeOptions, dict):
        raise SubEpiPredictError("SubEpiPredict:InvalidOdeOptions", "odeOptions must be an odeset structure.")
    if odeOptions.get("Events") is not None:
        raise SubEpiPredictError(
            "SubEpiPredict:ReservedEventsOption",
            "The simulator manages Events internally for subepidemic activation.",
        )
    if (
        odeOptions.get("Mass") is not None
        or odeOptions.get("Jacobian") is not None
        or str(odeOptions.get("Vectorized", "off")).lower() == "on"
    ):
        raise SubEpiPredictError(
            "SubEpiPredict:UnsupportedOdeOption",
            "Custom Mass/Jacobian and Vectorized='on' are not supported by this simulator.",
        )
    rtol = float(odeOptions.get("RelTol", 1e-8))
    atol = float(odeOptions.get("AbsTol", 1e-10))

    x, active, actTime, actCount, initialActive, nSegments = _simulate_core(
        t, float(I0), n, onset_fixed, float(onset_thr), flag1, rs, ps, as_, Ks, rtol, atol
    )
    totinc, patchIncidence = subepidemicIncidence(x, onset_fixed)
    activation = Struct(
        active=active, time=actTime, predecessorCount=actCount, initialActive=initialActive, nSegments=nSegments
    )
    return t, x, totinc, patchIncidence, activation


def subepidemicSimulateRefit(timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs1, ps1, as1, Ks1):
    """Simulate one bootstrap refit beyond the fitted window.

    Port of ``subepidemicSimulateRefit.m``. Returns
    ``(ok, x, totinc, patchIncidence, message)``. ``ok`` is False, and
    ``message`` says why, when the curve cannot be integrated to the last time
    or is not finite. Any other error (wrong sizes, invalid settings) is not
    caught.
    """
    try:
        _, x, totinc, patchIncidence, _ = simulateSubepidemic(
            timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs1, ps1, as1, Ks1
        )
    except SubEpiPredictError as simulationError:
        if simulationError.identifier not in INTEGRATION_FAILURES:
            raise
        return False, None, None, None, simulationError.message
    ok = x.shape[0] == np.size(timevect) and bool(np.all(np.isfinite(totinc)))
    message = "" if ok else "the curve is incomplete or not finite"
    return ok, x, totinc, patchIncidence, message
