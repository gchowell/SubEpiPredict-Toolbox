"""Settings of the fitting layer and the checks made before a fit.

Ports of ``subepidemicFitSettings.m`` (with ``options_fit_settings.m``),
``subepidemicCheckSettings.m``, ``subepidemicFitSeries.m``,
``subepidemicFitSignature.m`` and ``subepidemicSeed.m``.

MATLAB finds ``options_fit_settings.m`` on the path. In Python the overrides
are passed explicitly: every driver has a ``fit_settings`` argument that takes
a dict of overrides, a function that returns one (such as
:func:`options_fit_settings`), or the result of :func:`subepidemicFitSettings`.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError, subepidemicWarning
from .utils import Struct, get_rng, mat2str, movmean
from .utils import rng as _seed_rng

# fit revision of the MATLAB toolbox this port follows (subepidemicFitSignature.m)
# m2: bounds of r by kernel (rBoundsByKernel, rRateBounds); GGM, GLM and GRM fits are unchanged
FIT_REVISION = "m2"

_DEFAULTS = dict(
    logScaleR=True, logScaleK=True, logScaleAlpha=True,
    rLB=None, rUB=None, pLB=None, pUB=None, aLB=None, aUB=None, KLB=None, KUB=None,
    rUpperBoundFactor=10,
    rBoundsByKernel=True, rRateBounds=(1e-3, 10),
    countLikelihoodOnRawData=True,
    optimOptions=None,
    fitMaxAttempts=3,
    useParallel=False,
    polish=True, polishGradientTol=10, polishBudget=100,
    agreementTol=1e-2,
    continuationStart=True,
    bootstrapStarts=1,
    rngSeed=1,
    maxSkippedRefitFraction=0.05,
    baselineWindow=None,
    PIIscale=None,
    evidenceRatioLimit=100,
    firstThresholdStartFactor=3,
    notIdentifiedRatio=0.5,
)  # fmt: skip


def options_fit_settings():
    """Optional fit settings: the Python counterpart of ``options_fit_settings.m``.

    Returns a dict; each key replaces the default of the same name in
    :func:`subepidemicFitSettings`, every other setting keeps its default.
    Copy this function, uncomment what you need and pass it as
    ``fit_settings=`` to the drivers. As shipped it changes nothing. (The
    MATLAB file switches ``useParallel`` on; here it is off by default, see
    the README.)

    Examples of overrides::

        U['rLB'] = 1e-3                 # fixed floor for r
        U['rUB'] = 500                  # fixed ceiling for r
        U['rUpperBoundFactor'] = 0      # the legacy upper bound of r alone
        U['rBoundsByKernel'] = False    # count-based r bounds for every kernel, as before
        U['rRateBounds'] = (1e-4, 20)   # rate interval of r for logistic/Richards/Gompertz
        U['KUB'] = 5e6                  # carrying capacity at most 5e6
        U['pLB'] = 1; U['pUB'] = 1      # fix p of every sub-epidemic at 1
        U['KUB'] = [2e6, 5e6]           # one value per sub-epidemic
        U['useParallel'] = True         # local solves on a process pool
        U['optimOptions'] = {'FunctionTolerance': 1e-8}
        U['countLikelihoodOnRawData'] = False
        U['continuationStart'] = False
        U['bootstrapStarts'] = 3
        U['rngSeed'] = None             # new random numbers in every run
        U['maxSkippedRefitFraction'] = 0
        U['baselineWindow'] = 14
        U['evidenceRatioLimit'] = 1000
        U['firstThresholdStartFactor'] = 1
        U['notIdentifiedRatio'] = 0.8
    """
    U = Struct()
    return U


def _is_number(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def _empty(v):
    return v is None or (np.ndim(v) > 0 and np.size(v) == 0)


def subepidemicFitSettings(overrides=None):
    """Optimizer and search-space settings of the fitting layer.

    Port of ``subepidemicFitSettings.m``. ``overrides`` is a dict, a function
    returning a dict, or a settings Struct returned earlier by this function;
    each of its keys replaces the default of the same name. Unknown names are
    an error, so a misspelled setting cannot be ignored silently.

    SEARCH SPACE (:func:`subepidemicSearchSpace`)
      ``logScaleR``, ``logScaleK``, ``logScaleAlpha``  True: r, K and the
          negative-binomial alpha are searched on a log scale.
      ``rLB`` ... ``KUB``  None = the legacy bound. A scalar applies to every
          sub-epidemic; a vector gives one value per sub-epidemic. A lower
          bound equal to its upper bound fixes the parameter.
      ``rUpperBoundFactor``  10: the default upper bound of r is the larger of
          50*max(first two observations) and 10 times the largest count.
          0 gives the legacy bound alone.
      ``rBoundsByKernel``  True: the kernels with p = 1 (logistic, Richards,
          Gompertz), whose r is a rate per time step, take the bounds of r
          from ``rRateBounds`` instead of the counts: both ends for logistic
          and Gompertz, the lower end for Richards (r trades off against a).
          False gives the count-based bounds for every kernel, as before.
      ``rRateBounds``  (1e-3, 10): [lower, upper] rate per time step,
          0 < lower <= upper. ``rLB`` and ``rUB`` still take precedence.

    FITTED SERIES (:func:`subepidemicFitSeries`)
      ``countLikelihoodOnRawData``  True: Poisson and negative-binomial fits
          use the raw counts whatever ``smoothfactor1`` is.

    OPTIMIZER (:func:`fitSubepidemicModel`)
      ``optimOptions``  dict of overrides with the fmincon names
          ``FunctionTolerance``, ``FiniteDifferenceStepSize``,
          ``FiniteDifferenceType``, ``MaxIterations``,
          ``MaxFunctionEvaluations`` (``StepTolerance`` and
          ``OptimalityTolerance`` are accepted and recorded).
      ``fitMaxAttempts``  3 multi-start passes before a candidate fails.
      ``useParallel``  False. True runs the local solves on a process pool.
      ``polish``, ``polishGradientTol`` (10), ``polishBudget`` (100)
          Nelder-Mead pass, run only when the stationarity check fails.
      ``agreementTol``  0.01 log-likelihood units.

    SCREEN AND BOOTSTRAP
      ``continuationStart``  True; ``bootstrapStarts``  1;
      ``firstThresholdStartFactor``  3.

    REPEATABLE RUNS
      ``rngSeed``  1; None leaves the generator alone.

    REPORTS
      ``maxSkippedRefitFraction`` 0.05, ``baselineWindow`` None,
      ``PIIscale`` None, ``evidenceRatioLimit`` 100,
      ``notIdentifiedRatio`` 0.5.
    """
    S = Struct(_DEFAULTS)
    S.optimOptions = Struct()

    U = overrides
    if callable(U):
        U = U()
    if U is not None:
        if not isinstance(U, dict):
            raise SubEpiPredictError(
                "SubEpiPredict:FitSettings", "options_fit_settings must return one struct (dict) of settings."
            )
        for name, value in U.items():
            if name not in S:
                raise SubEpiPredictError(
                    "SubEpiPredict:FitSettings",
                    "Unknown fit setting '%s'. See subepidemicFitSettings for the list." % name,
                )
            S[name] = value
    if S.optimOptions is None:
        S.optimOptions = Struct()

    # ---- validation -------------------------------------------------------------
    for name in ("logScaleR", "logScaleK", "logScaleAlpha", "countLikelihoodOnRawData", "useParallel", "polish",
                 "continuationStart", "rBoundsByKernel"):  # fmt: skip
        v = S[name]
        if not ((isinstance(v, (bool, np.bool_)) or _is_number(v)) and (v == 0 or v == 1)):
            raise SubEpiPredictError("SubEpiPredict:FitSettings", "Setting %s must be true or false." % name)
        S[name] = bool(v)
    for name in ("polishGradientTol", "polishBudget", "agreementTol"):
        v = S[name]
        if not (_is_number(v) and math.isfinite(v) and v > 0):
            raise SubEpiPredictError("SubEpiPredict:FitSettings", "Setting %s must be a positive number." % name)
    v = S.rUpperBoundFactor
    if not (_is_number(v) and math.isfinite(v) and v >= 0):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting rUpperBoundFactor must be a number of at least 0."
        )
    try:
        rb = np.asarray(S.rRateBounds, dtype=float).ravel()
        ok = rb.size == 2 and bool(np.all(np.isfinite(rb))) and 0 < rb[0] <= rb[1]
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting rRateBounds must be [lower upper] with 0 < lower <= upper."
        )
    S.rRateBounds = (float(rb[0]), float(rb[1]))
    for name, minimum in (("fitMaxAttempts", 1), ("bootstrapStarts", 0)):
        v = S[name]
        if not (_is_number(v) and math.isfinite(v) and v == int(v) and v >= minimum):
            raise SubEpiPredictError(
                "SubEpiPredict:FitSettings", "Setting %s must be an integer of at least %d." % (name, minimum)
            )
        S[name] = int(v)
    for name in ("rLB", "rUB", "pLB", "pUB", "aLB", "aUB", "KLB", "KUB"):
        v = S[name]
        if _empty(v):
            S[name] = None
            continue
        try:
            a = np.asarray(v, dtype=float)
            ok = a.ndim <= 1 and bool(np.all(np.isfinite(a)))
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise SubEpiPredictError(
                "SubEpiPredict:FitSettings",
                "Setting %s must be empty, a finite scalar, or a finite vector." % name,
            )
    v = S.rngSeed
    if _empty(v):
        S.rngSeed = None
    elif not (_is_number(v) and math.isfinite(v) and v >= 0 and v == int(v) and v < 2**32):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting rngSeed must be empty or a nonnegative integer below 2^32."
        )
    else:
        S.rngSeed = int(v)
    v = S.baselineWindow
    if _empty(v):
        S.baselineWindow = None
    elif not (_is_number(v) and math.isfinite(v) and v >= 1 and v == int(v)):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting baselineWindow must be empty or a positive integer."
        )
    v = S.maxSkippedRefitFraction
    if not (_is_number(v) and 0 <= v <= 1):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting maxSkippedRefitFraction must be a number from 0 to 1."
        )
    v = S.firstThresholdStartFactor
    if not (_is_number(v) and math.isfinite(v) and v >= 1):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings",
            "Setting firstThresholdStartFactor must be a number of at least 1 (1: no extra start points).",
        )
    v = S.notIdentifiedRatio
    if not (_is_number(v) and not math.isnan(v) and v > 0):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings",
            "Setting notIdentifiedRatio must be a positive number (Inf turns the warning off).",
        )
    v = S.evidenceRatioLimit
    if not (_is_number(v) and not math.isnan(v) and v > 1):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings",
            "Setting evidenceRatioLimit must be a number above 1 (Inf turns the warning off).",
        )
    v = S.PIIscale
    if _empty(v):
        S.PIIscale = None
    else:
        try:
            a = np.asarray(v, dtype=float).ravel()
            ok = a.size == 6 and bool(np.all(np.isnan(a) | (a == 0) | (a == 1)))
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise SubEpiPredictError(
                "SubEpiPredict:FitSettings",
                "Setting PIIscale must be empty or 6 entries [r p a K alpha d], each 1 (log scale), "
                "0 (range fraction) or NaN (automatic).",
            )
    if not isinstance(S.optimOptions, dict):
        raise SubEpiPredictError(
            "SubEpiPredict:FitSettings", "Setting optimOptions must be a struct (dict) of optimizer options."
        )
    return S


def resolve_fit_settings(fit_settings=None):
    """The settings Struct for a ``fit_settings`` argument (None, dict, function or Struct)."""
    if isinstance(fit_settings, Struct) and set(_DEFAULTS).issubset(fit_settings.keys()):
        return fit_settings
    return subepidemicFitSettings(fit_settings)


# <============================================================================>
# <============================ Check of the options ==========================>
# <============================================================================>


def _show(v):
    if isinstance(v, str):
        return "'" + v + "'"
    try:
        a = np.asarray(v, dtype=float)
    except (TypeError, ValueError):
        return "a " + type(v).__name__
    if a.size == 0:
        return "[]"
    if a.size <= 6:
        return mat2str(a, 6)
    return "[%d values]" % a.size


def _join(problems):
    if len(problems) == 1:
        return "The settings cannot run. " + problems[0]
    s = "The settings cannot run (%d problems):" % len(problems)
    for k, problem in enumerate(problems, start=1):
        s = "%s\n  %d. %s" % (s, k, problem)
    return s


def subepidemicCheckSettings(S):
    """Check the settings of a fit before any fitting starts.

    Port of ``subepidemicCheckSettings.m``. Returns silently when the settings
    can run and otherwise raises ``SubEpiPredict:InvalidSettings`` with a
    message that lists EVERY problem found.

    ``S`` is a dict (or an :class:`~subepipredict.Options`) with the
    fields ``method1``, ``dist1``, ``flag1``, ``npatches_fixed``,
    ``onset_fixed``, ``smoothfactor1``, ``calibrationperiod1``,
    ``numstartpoints``, ``topmodelsx``, ``B``, ``DT``, ``cumulative1`` and
    ``outbreakx``.

    Estimation method and error structure must match: method 0 with dist 0, 1
    or 2; method 1 with dist 1; methods 3, 4, 5 with the same dist; method 6
    with dist 6 or 1.
    """
    if hasattr(S, "as_dict"):
        S = S.as_dict()
    if not isinstance(S, dict):
        raise SubEpiPredictError("SubEpiPredict:InvalidSettings", "The settings must be passed as one struct.")
    names = ["method1", "dist1", "flag1", "npatches_fixed", "onset_fixed", "smoothfactor1", "calibrationperiod1",
             "numstartpoints", "topmodelsx", "B", "DT", "cumulative1", "outbreakx"]  # fmt: skip
    problems = ["%s is missing." % name for name in names if name not in S]
    if problems:
        raise SubEpiPredictError("SubEpiPredict:InvalidSettings", _join(problems))

    def isInt(v):
        return (
            (_is_number(v) or isinstance(v, (bool, np.bool_))) and math.isfinite(v) and v == int(v)
        )

    # ---- estimation method and error structure ----
    methodOK = isInt(S["method1"]) and S["method1"] in (0, 1, 3, 4, 5, 6)
    if not methodOK:
        if isInt(S["method1"]) and S["method1"] == 2:
            problems.append("method1 = 2 (Pearson chi-squared) is not implemented. Use 0, 1, 3, 4, 5 or 6.")
        else:
            problems.append("method1 = %s is not valid. Use 0, 1, 3, 4, 5 or 6." % _show(S["method1"]))
    distOK = isInt(S["dist1"]) and S["dist1"] in (0, 1, 2, 3, 4, 5, 6)
    if not distOK:
        problems.append("dist1 = %s is not valid. Use an integer from 0 to 6." % _show(S["dist1"]))
    if methodOK and distOK:
        m = int(S["method1"])
        if m == 0:
            allowed = (0, 1, 2)
            hint = (
                "Least squares (method1 = 0) needs dist1 = 0 (normal), 1 (Poisson) or 2 "
                "(negative binomial with VAR = factor*mean). The options do not set dist1 for method1 = 0."
            )
        elif m == 1:
            allowed = (1,)
            hint = "The Poisson likelihood (method1 = 1) needs dist1 = 1."
        elif m in (3, 4, 5):
            allowed = (m,)
            hint = "The negative binomial likelihood method1 = %d needs dist1 = %d." % (m, m)
        else:
            allowed = (6, 1)
            hint = "Least absolute deviations (method1 = 6) needs dist1 = 6 (Laplace) or 1 (Poisson)."
        if S["dist1"] not in allowed:
            problems.append("method1 = %d with dist1 = %d cannot run. %s" % (m, S["dist1"], hint))
    # ---- growth model ----
    if not (isInt(S["flag1"]) and 0 <= S["flag1"] <= 5):
        problems.append(
            "flag1 = %s is not valid. Use one integer from 0 to 5 for all sub-epidemics; a vector is not supported."
            % _show(S["flag1"])
        )
    # ---- sizes and switches ----
    if not (isInt(S["npatches_fixed"]) and S["npatches_fixed"] >= 1):
        problems.append("npatches_fixed = %s must be a positive integer." % _show(S["npatches_fixed"]))
    if not (isInt(S["onset_fixed"]) and S["onset_fixed"] in (0, 1)):
        problems.append("onset_fixed = %s must be 0 or 1." % _show(S["onset_fixed"]))
    if not (isInt(S["smoothfactor1"]) and S["smoothfactor1"] >= 1):
        problems.append(
            "smoothfactor1 = %s must be a positive integer (1 = no smoothing)." % _show(S["smoothfactor1"])
        )
    if not (isInt(S["calibrationperiod1"]) and S["calibrationperiod1"] >= 2):
        problems.append("calibrationperiod1 = %s must be an integer of at least 2." % _show(S["calibrationperiod1"]))
    if not (isInt(S["numstartpoints"]) and S["numstartpoints"] >= 1):
        problems.append("numstartpoints = %s must be a positive integer." % _show(S["numstartpoints"]))
    if not (isInt(S["topmodelsx"]) and S["topmodelsx"] >= 1):
        problems.append("topmodelsx = %s must be a positive integer." % _show(S["topmodelsx"]))
    if not (isInt(S["B"]) and S["B"] >= 1):
        problems.append("B = %s (bootstrap realizations) must be a positive integer." % _show(S["B"]))
    if not (isInt(S["DT"]) and S["DT"] in (1, 7, 365)):
        problems.append("DT = %s must be 1 (daily), 7 (weekly) or 365 (yearly)." % _show(S["DT"]))
    if not (isInt(S["cumulative1"]) and S["cumulative1"] in (0, 1)):
        problems.append(
            "cumulative1 = %s must be 0 (incidence file) or 1 (cumulative file)." % _show(S["cumulative1"])
        )
    v = S["outbreakx"]
    try:
        a = np.asarray(v, dtype=float)
        ok = a.size > 0 and a.ndim <= 1 and bool(np.all(np.isfinite(a)) and np.all(a == np.trunc(a)) and np.all(a >= 1))
        ok = ok and not isinstance(v, str)
    except (TypeError, ValueError):
        ok = False
    if not ok:
        problems.append("outbreakx = %s must hold positive integers (columns of the data file)." % _show(v))
    if problems:
        raise SubEpiPredictError("SubEpiPredict:InvalidSettings", _join(problems))


# <============================================================================>
# <================================ Fitted series =============================>
# <============================================================================>


def subepidemicFitSeries(data, smoothfactor1, method1, settings=None):
    """The series the objective is fitted to.

    Port of ``subepidemicFitSeries.m``. The screen and the bootstrap both call
    this function, so they cannot disagree about what was fitted.

    Least squares (method 0) and absolute deviations (method 6): the moving
    average of ``data`` with span ``smoothfactor1``. Count likelihoods
    (methods 1, 3, 4, 5): the raw counts, whatever ``smoothfactor1`` is,
    unless ``settings.countLikelihoodOnRawData`` is False.

    Count likelihoods need nonnegative data: a negative value is an error.
    Non-integer raw values are accepted with a warning.

    Returns ``(ydata, info)``; ``info`` has ``usedRawData``,
    ``smoothingSpan``, ``smoothedCountLikelihood`` and ``nNonInteger``.
    """
    settings = resolve_fit_settings(settings)
    try:
        data = np.asarray(data, dtype=float)
        ok = data.size > 0 and (data.ndim <= 1 or (data.ndim == 2 and 1 in data.shape))
        ok = ok and bool(np.all(np.isfinite(data)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError("SubEpiPredict:InvalidData", "data must be a nonempty finite real vector.")
    if not (_is_number(method1) and method1 in (0, 1, 2, 3, 4, 5, 6)):
        raise SubEpiPredictError("SubEpiPredict:InvalidMethod", "method1 must be an integer from 0 to 6.")
    if _empty(smoothfactor1):
        smoothfactor1 = 1
    if not (_is_number(smoothfactor1) and math.isfinite(smoothfactor1) and smoothfactor1 >= 1):
        raise SubEpiPredictError("SubEpiPredict:InvalidSmoothing", "smoothfactor1 must be a number of at least 1.")
    data = data.ravel().copy()
    isCount = method1 in (1, 3, 4, 5)
    useRaw = isCount and settings.countLikelihoodOnRawData
    info = Struct(usedRawData=False, smoothingSpan=1, smoothedCountLikelihood=False, nNonInteger=0)
    if useRaw or smoothfactor1 <= 1:
        ydata = data
        info.usedRawData = True
    else:
        ydata = movmean(data, smoothfactor1)
        info.smoothingSpan = smoothfactor1
        info.smoothedCountLikelihood = isCount
    if isCount and info.usedRawData:
        nNegative = int(np.sum(data < 0))
        if nNegative > 0:
            raise SubEpiPredictError(
                "SubEpiPredict:NegativeCounts",
                "The calibration series has %d negative value(s), and method1 = %d fits a count likelihood "
                "to the raw series. Correct the data, or use least squares (method1 = 0), or set "
                "countLikelihoodOnRawData = False in the fit settings to fit the moving average as before."
                % (nNegative, method1),
            )
        info.nNonInteger = int(np.sum(data != np.round(data)))
        if info.nNonInteger > 0:
            subepidemicWarning(
                "SubEpiPredict:NonIntegerCounts",
                "The calibration series has %d non-integer value(s) and method1 = %d fits a count likelihood. "
                "The likelihood keeps its convention for fractional data; integer counts are what it describes.",
                info.nNonInteger, method1,
            )  # fmt: skip
    return ydata, info


# <============================================================================>
# <============================== Fit signature ===============================>
# <============================================================================>

_REPORT_ONLY = ("PIIscale", "evidenceRatioLimit", "notIdentifiedRatio", "baselineWindow", "maxSkippedRefitFraction")


def _same(a, b):
    """isequaln for the values a setting can hold."""
    if a is None or b is None:
        return _empty(a) and _empty(b)
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)) or set(a) != set(b):
            return False
        return all(_same(a[k], b[k]) for k in a)
    if isinstance(a, str) or isinstance(b, str):
        return str(a) == str(b)
    try:
        aa = np.asarray(a, dtype=float).ravel()
        bb = np.asarray(b, dtype=float).ravel()
    except (TypeError, ValueError):
        return a == b
    return aa.size == bb.size and bool(np.all((aa == bb) | (np.isnan(aa) & np.isnan(bb))))


def subepidemicFitSignature(saved=None, settings=None, compare=False):
    """What made a fit: the fit revision and the settings that shape it.

    Port of ``subepidemicFitSignature.m``.

    ``subepidemicFitSignature(settings=...)`` returns the signature of a fit
    made now: a Struct with ``fitRevision`` and ``settings`` (the fit settings
    without those that only change reports).

    ``subepidemicFitSignature(saved, settings=..., compare=True)`` returns
    ``(same, differences)``: whether a fit made now would have the signature
    ``saved``, and the names of what differs. ``saved = None`` stands for a
    fit stored without a signature; ``same`` is then False.
    """
    settings = resolve_fit_settings(settings)
    shaped = Struct({k: v for k, v in settings.items() if k not in _REPORT_ONLY})
    signature = Struct(fitRevision=FIT_REVISION, settings=shaped)
    if not compare:
        return signature

    differences = []
    if not (isinstance(saved, dict) and "fitRevision" in saved and isinstance(saved.get("settings"), dict)):
        differences = ["no signature (fit made before revision m1)"]
    else:
        if str(saved["fitRevision"]) != signature.fitRevision:
            differences.append("fit revision (%s then, %s now)" % (saved["fitRevision"], signature.fitRevision))
        for name in sorted(set(saved["settings"]) | set(signature.settings)):
            if (
                name not in saved["settings"]
                or name not in signature.settings
                or not _same(saved["settings"][name], signature.settings[name])
            ):
                differences.append(name)
    return len(differences) == 0, differences


# <============================================================================>
# <================================ Random seed ===============================>
# <============================================================================>


def subepidemicSeed(label=None, settings=None, verbose=True):
    """Apply the random seed setting and record it.

    Port of ``subepidemicSeed.m``. Called at the start of the fit driver (once
    per area) and of the forecast, fit-plot and Rt functions. With
    ``settings.rngSeed`` an integer, the toolbox stream is reset with that
    seed, so the same call on the same files gives the same numbers again;
    with None the stream is left as it is.

    Returns a Struct with ``seed`` (None when no seed was set) and ``state``
    (the state of the generator at this point, as returned by
    ``Generator.bit_generator.state``).
    """
    if not label:
        label = "SubEpiPredict"
    settings = resolve_fit_settings(settings)
    info = Struct(seed=settings.rngSeed, state=None)
    if info.seed is not None:
        _seed_rng(info.seed)
        if verbose:
            print("%s: random seed %d (setting rngSeed)." % (label, info.seed))
    info.state = get_rng().bit_generator.state
    return info
