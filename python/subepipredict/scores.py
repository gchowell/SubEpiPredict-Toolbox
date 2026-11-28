"""Forecast scores, quantiles, ensemble weights and the flat baseline.

Ports of ``computeQuantiles.m``, ``computeWIS.m``,
``computeforecastperformance.m``, ``subepidemicWeightedQuantile.m``,
``subepidemicBaselineForecast.m``, ``subepidemicEnsembleWeights.m`` and
``subepidemicEvidenceRatios.m``.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError, subepidemicWarning
from .io import ALPHAQUANTILES
from .utils import Struct, mat2str, quantile

# central prediction intervals of the weighted interval score (COVID-19 Forecast Hub)
WIS_ALPHAS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def _is_real_scalar(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def _clip0(q):
    """``(q >= 0) .* q``: negative values become 0, NaN stays NaN."""
    q = np.asarray(q, dtype=float)
    with np.errstate(invalid="ignore"):
        return (q >= 0) * q + 0.0  # + 0.0 turns -0.0 into 0.0


# <============================================================================>
# <============================ Weighted quantiles ============================>
# <============================================================================>


def subepidemicWeightedQuantile(X, p, w=None):
    """Quantiles of weighted curves at each time point.

    Port of ``subepidemicWeightedQuantile.m``.

    ``X``  T-by-N, one row per time point and one column per curve
    ``p``  probabilities, each from 0 to 1
    ``w``  N nonnegative weights, one per curve; they need not sum to 1.
           None = every curve weighs the same.

    Returns ``q``, T-by-len(p); ``q[i, k]`` is the ``p[k]`` quantile of row
    ``i`` of ``X``.

    Definition. At one time point, sort the values and carry their weights
    along, scaled to sum to 1. The value ``x(j)`` is the quantile of
    probability ``P(j) = w(1) + ... + w(j-1) + w(j)/2``; quantiles between two
    values are linearly interpolated, and probabilities below ``P(1)`` or
    above ``P(n)`` give ``x(1)`` or ``x(n)``. With equal weights this is the
    definition of MATLAB's ``quantile``. Equal values share their weights
    equally, so the result does not depend on the order of the curves. Curves
    with weight 0 and NaN values are left out; a row with no value left gives
    NaN.
    """
    try:
        X = np.asarray(X, dtype=float)
    except (TypeError, ValueError):
        X = None
    if X is not None and X.ndim == 1:
        X = X.reshape(1, -1)
    if X is None or X.ndim != 2:
        raise SubEpiPredictError(
            "SubEpiPredict:WeightedQuantile",
            "X must be a real matrix, one row per time point and one column per curve.",
        )
    try:
        p = np.asarray(p, dtype=float).ravel()
        ok = p.size > 0 and bool(np.all((p >= 0) & (p <= 1)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError("SubEpiPredict:WeightedQuantile", "The probabilities must be numbers from 0 to 1.")
    T, N = X.shape
    if w is None or np.size(w) == 0:
        w = np.ones(N)
    try:
        w = np.asarray(w, dtype=float).ravel()
        ok = w.size == N and bool(np.all(np.isfinite(w)) and np.all(w >= 0))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:WeightedQuantile",
            "Give one finite, nonnegative weight per curve (%d curves, %d weights)." % (N, np.size(w)),
        )
    q = np.full((T, p.size), np.nan)
    used = w > 0
    for i in range(T):
        x = X[i, :]
        ok = used & ~np.isnan(x)
        if not np.any(ok):
            continue
        order = np.argsort(x[ok], kind="stable")
        xs = x[ok][order]
        ws = w[ok][order]
        ws = ws / ws.sum()
        n = xs.size
        # equal values share their weights equally
        group = np.concatenate(([0], np.cumsum(np.diff(xs) != 0)))
        if group[-1] < n - 1:
            groupWeight = np.bincount(group, weights=ws)
            groupSize = np.bincount(group)
            ws = groupWeight[group] / groupSize[group]
        # probability of each sorted value: the weight below it plus half its own
        P = np.cumsum(ws) - ws / 2
        j = np.searchsorted(P, p, side="right")  # number of values at or below p
        out = np.empty(p.size)
        low = j == 0
        high = j == n
        mid = ~(low | high)
        out[low] = xs[0]
        out[high] = xs[n - 1]
        if np.any(mid):
            jm = j[mid]
            # P(j) <= p < P(j+1), so the denominator is positive
            out[mid] = xs[jm - 1] + (p[mid] - P[jm - 1]) / (P[jm] - P[jm - 1]) * (xs[jm] - xs[jm - 1])
        q[i, :] = out
    return q


def _rowQuantile(curves, p, weights=None):
    """Quantiles ``p`` of the curves (columns) at each time point (row): rows-by-len(p)."""
    curves = np.asarray(curves, dtype=float)
    if curves.ndim == 1:
        curves = curves.reshape(1, -1)
    p = np.atleast_1d(np.asarray(p, dtype=float))
    if weights is None or np.size(weights) == 0:
        if curves.shape[0] == 0:
            return np.zeros((0, p.size))
        return np.asarray(quantile(curves, p, axis=1)).T
    return subepidemicWeightedQuantile(curves, p, weights)


def _rowMedian(curves, weights=None):
    return _rowQuantile(curves, [0.5], weights)[:, 0]


# <============================================================================>
# <================================ Quantiles =================================>
# <============================================================================>


def computeQuantiles(data1, curvesforecasts2, forecastingperiod, weights2=None):
    """The 23 predictive quantiles of the calibration and the forecast period.

    Port of ``computeQuantiles.m``. ``data1`` has one row per calibration time
    point; ``curvesforecasts2`` one row per time point (calibration followed
    by forecast) and one column per curve; ``weights2`` (optional) one weight
    per curve.

    Returns ``(quantilesc1, quantilesf1, alphaquantiles)``; negative quantiles
    are set to 0.
    """
    curves = np.asarray(curvesforecasts2, dtype=float)
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)
    quantilesc1 = _clip0(_rowQuantile(curves[:tf2, :], ALPHAQUANTILES, weights2))
    quantilesf1 = np.zeros((forecastingperiod, len(ALPHAQUANTILES)))
    if forecastingperiod > 0:
        quantilesf1 = _clip0(_rowQuantile(curves[tf2 : tf2 + forecastingperiod, :], ALPHAQUANTILES, weights2))
    return quantilesc1, quantilesf1, np.array(ALPHAQUANTILES)


# <============================================================================>
# <========================= Weighted interval score ==========================>
# <============================================================================>


def _wis_rows(y, curves, weights):
    """WIS of each observation ``y[i]`` against row ``i`` of ``curves``."""
    alphas = np.array(WIS_ALPHAS)
    K = alphas.size
    w0 = 1 / 2
    probs = np.concatenate([alphas / 2, 1 - alphas / 2, [0.5]])
    Q = _rowQuantile(curves, probs, weights)
    Lt = _clip0(Q[:, :K])
    Ut = _clip0(Q[:, K : 2 * K])
    m = Q[:, 2 * K]
    y = np.asarray(y, dtype=float).reshape(-1, 1)
    IS = (Ut - Lt) + (2 / alphas) * (Lt - y) * (y < Lt) + (2 / alphas) * (y - Ut) * (y > Ut)
    sum1 = np.sum((alphas / 2) * IS, axis=1)
    return (1 / (K + 1 / 2)) * (w0 * np.abs(y[:, 0] - m) + sum1)


def computeWIS(data1, datalatest, curvesforecasts2, forecastingperiod, weights2=None):
    """Weighted interval score of the calibration period and of every horizon.

    Port of ``computeWIS.m``. The score uses the median and the 11 central
    intervals of the COVID-19 Forecast Hub (alpha = 0.02, 0.05, 0.1, ..., 0.9).

    ``data1``             calibration data (rows = time points)
    ``datalatest``        ``[time, observation]`` for the calibration period
                          followed by the forecast period
    ``curvesforecasts2``  predictive curves, rows = time points, columns = curves
    ``forecastingperiod`` number of forecast time points to score (0: none)
    ``weights2``          optional weight per curve (ensembles)

    Returns ``(WISC, WISFS, quantilesc1, quantilesf1)``: ``WISC`` the mean WIS
    of the calibration period; ``WISFS`` one row ``[h, mean WIS over the first
    h forecast time points]`` per horizon; the 23 quantiles of both periods.
    """
    curves = np.asarray(curvesforecasts2, dtype=float)
    datalatest = np.asarray(datalatest, dtype=float)
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)

    # calibration period
    datac = datalatest[:tf2, 1]
    WISC = float(np.mean(_wis_rows(datac, curves[:tf2, :], weights2)))
    quantilesc1 = _clip0(_rowQuantile(curves[:tf2, :], ALPHAQUANTILES, weights2))

    # forecasting period
    WISFS = np.zeros((0, 2))
    quantilesf1 = np.zeros((forecastingperiod, len(ALPHAQUANTILES)))
    if forecastingperiod > 0:
        rows = slice(tf2, tf2 + forecastingperiod)
        dataf = datalatest[rows, 1]
        WISF = _wis_rows(dataf, curves[rows, :], weights2)
        horizons = np.arange(1, forecastingperiod + 1)
        WISFS = np.column_stack([horizons, np.cumsum(WISF) / horizons])
        quantilesf1 = _clip0(_rowQuantile(curves[rows, :], ALPHAQUANTILES, weights2))
    return WISC, WISFS, quantilesc1, quantilesf1


# <============================================================================>
# <============================ Performance metrics ===========================>
# <============================================================================>


def computeforecastperformance(
    data1, datalatest, curvesforecasts1, curvesforecasts2, forecastingperiod, weights1=None, weights2=None
):
    """Error, coverage and interval score of the calibration and forecast periods.

    Port of ``computeforecastperformance.m``.

    ``curvesforecasts1``  model trajectories (their median is the point forecast)
    ``curvesforecasts2``  predictive curves (trajectory plus noise), for the
                          95% prediction interval
    ``weights1``, ``weights2``  optional weight per column (ensembles)

    Returns ``(RMSECS, MSECS, MAECS, PICS, MISCS, RMSEFS, MSEFS, MAEFS, PIFS,
    MISFS)``. The first five are the numbers of the calibration period (root
    mean squared error, mean squared error, mean absolute error, coverage of
    the 95% prediction interval in percent, mean interval score). The last
    five have one row ``[h, value over the first h forecast time points]`` per
    horizon (empty when ``forecastingperiod`` is 0).
    """
    curves1 = np.asarray(curvesforecasts1, dtype=float)
    curves2 = np.asarray(curvesforecasts2, dtype=float)
    if curves1.ndim == 1:
        curves1 = curves1.reshape(-1, 1)
    if curves2.ndim == 1:
        curves2 = curves2.reshape(-1, 1)
    datalatest = np.asarray(datalatest, dtype=float)
    tf2 = np.asarray(data1).shape[0]
    forecastingperiod = int(forecastingperiod)

    def scores(y, c1, c2):
        ymean = _rowMedian(c1, weights1)
        Q = _rowQuantile(c2, [0.025, 0.975], weights2)
        Lt = _clip0(Q[:, 0])
        Ut = _clip0(Q[:, 1])
        err = y - ymean
        inside = (y >= Lt) & (y <= Ut)
        IS = (Ut - Lt) + (2 / 0.05) * (Lt - y) * (y < Lt) + (2 / 0.05) * (y - Ut) * (y > Ut)
        return err, inside, IS

    # calibration period
    datac = datalatest[:tf2, 1]
    err, inside, IS = scores(datac, curves1[:tf2, :], curves2[:tf2, :])
    MSEc = float(np.mean(err**2))
    RMSEc = math.sqrt(MSEc)
    MAEc = float(np.mean(np.abs(err)))
    PIc = 100 * float(np.sum(inside)) / datac.size
    MISc = float(np.mean(IS))

    # forecasting period
    empty = np.zeros((0, 2))
    RMSEFS, MSEFS, MAEFS, PIFS, MISFS = empty, empty, empty, empty, empty
    if forecastingperiod > 0:
        rows = slice(tf2, tf2 + forecastingperiod)
        dataf = datalatest[rows, 1]
        err, inside, IS = scores(dataf, curves1[rows, :], curves2[rows, :])
        h = np.arange(1, forecastingperiod + 1)
        mse = np.cumsum(err**2) / h
        MSEFS = np.column_stack([h, mse])
        RMSEFS = np.column_stack([h, np.sqrt(mse)])
        MAEFS = np.column_stack([h, np.cumsum(np.abs(err)) / h])
        PIFS = np.column_stack([h, 100 * np.cumsum(inside) / h])
        MISFS = np.column_stack([h, np.cumsum(IS) / h])
    return RMSEc, MSEc, MAEc, PIc, MISc, RMSEFS, MSEFS, MAEFS, PIFS, MISFS


# <============================================================================>
# <=============================== Flat baseline ==============================>
# <============================================================================>


def _sortedQuantile(x, p):
    """Quantiles of the sorted vector x by the rule of MATLAB's quantile."""
    n = x.size
    position = np.minimum(np.maximum(p * n + 0.5, 1), n)
    lower = np.floor(position).astype(int)
    upper = np.minimum(lower + 1, n)
    weight = position - lower
    return x[lower - 1] * (1 - weight) + x[upper - 1] * weight


def subepidemicBaselineForecast(y, horizon, window, nSamples=None):
    """Flat baseline forecast with empirical uncertainty.

    Port of ``subepidemicBaselineForecast.m``. The reference a model forecast
    should beat.

    Point forecast: the mean of the last ``window`` observations, for every
    horizon. Uncertainty: for horizon h the same rule is applied at every
    earlier time t of the series and its error
    ``e = y(t+h) - mean(y(t-window+1 : t))`` recorded; the forecast
    distribution at horizon h is the point forecast plus these errors and
    their negatives (both signs keep the median at the point forecast).
    Values below 0 are set to 0. If fewer than 5 errors exist for a horizon,
    the errors of the longest shorter horizon that has 5 are used; if none
    has, the uncertainty is not defined (NaN).

    Returns a Struct: ``level`` (the point forecast), ``point`` (one value per
    horizon), ``samples`` (horizon-by-nSamples, default 200: evenly spaced
    quantiles of the forecast distribution of each horizon; no random
    numbers), ``nErrors``, ``borrowed`` and ``window``.
    """
    if nSamples is None:
        nSamples = 200
    try:
        y = np.asarray(y, dtype=float)
        ok = y.size > 0 and (y.ndim <= 1 or (y.ndim == 2 and 1 in y.shape)) and bool(np.all(np.isfinite(y)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError("SubEpiPredict:Baseline", "The series must be a non-empty vector of finite numbers.")
    if not (_is_real_scalar(horizon) and math.isfinite(horizon) and horizon >= 1 and horizon == int(horizon)):
        raise SubEpiPredictError("SubEpiPredict:Baseline", "The horizon must be a positive integer.")
    if not (_is_real_scalar(window) and math.isfinite(window) and window >= 1 and window == int(window)):
        raise SubEpiPredictError("SubEpiPredict:Baseline", "The window must be a positive integer.")
    y = y.ravel()
    n = y.size
    horizon = int(horizon)
    nSamples = int(nSamples)
    window = min(int(window), n)
    B = Struct(
        level=float(np.mean(y[n - window :])), point=None, samples=np.full((horizon, nSamples), np.nan),
        nErrors=np.zeros(horizon, dtype=int), borrowed=np.zeros(horizon, dtype=bool), window=window,
    )  # fmt: skip
    B.point = np.full(horizon, B.level)

    # rolling mean of the last `window` values, available from t = window on
    c = np.concatenate(([0.0], np.cumsum(y)))
    level = (c[window : n + 1] - c[0 : n - window + 1]) / window  # level[k] is the mean ending at t = window+k (1-based)
    minErrors = 5
    probabilities = (np.arange(1, nSamples + 1) - 0.5) / nSamples
    lastGood = None
    for h in range(1, horizon + 1):
        nOrigins = n - h - window + 1  # origins t = window, ..., n-h
        if nOrigins >= 1:
            e = y[window + h - 1 : n] - level[:nOrigins]
        else:
            e = np.zeros(0)
        B.nErrors[h - 1] = e.size
        if e.size >= minErrors:
            lastGood = e
        elif lastGood is not None:
            e = lastGood
            B.borrowed[h - 1] = True
            B.nErrors[h - 1] = e.size
        else:
            continue  # uncertainty not defined at this horizon
        both = np.sort(np.concatenate([e, -e]))
        B.samples[h - 1, :] = np.maximum(B.level + _sortedQuantile(both, probabilities), 0)
    return B


# <============================================================================>
# <============================== Ensemble weights ============================>
# <============================================================================>


def _inverseScore(score, k, name, code):
    try:
        score = np.asarray(score, dtype=float).ravel()
        ok = score.size == k
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleWeights",
            "weight_type1 = %d needs one %s value per member model (%d)." % (code, name, k),
        )
    if not (np.all(np.isfinite(score)) and np.all(score > 0)):
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleWeights",
            "weight_type1 = %d needs finite, positive %s values; received %s." % (code, name, mat2str(score, 6)),
        )
    return (1 / score) / np.sum(1 / score)


def subepidemicEnsembleWeights(weight_type1, AICc, WISC=None, WISF=None):
    """Weights of the member models of an ensemble.

    Port of ``subepidemicEnsembleWeights.m``. ``AICc`` holds the AICc of the
    member models; ``WISC`` and ``WISF`` their weighted interval scores
    (calibration; earlier forecasts), read only for codes 2 and 3. The result
    holds nonnegative weights that sum to 1.

    ``weight_type1``: -1 equal weights; 1 Akaike weights,
    ``exp(-(AICc - min(AICc))/2)`` normalized; 0 the same with a warning (it
    used to be ``1/AICc``, which is not a weight); 2 proportional to
    ``1/WISC``; 3 proportional to ``1/WISF``.
    """
    if not (_is_real_scalar(weight_type1) and weight_type1 in (-1, 0, 1, 2, 3)):
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleWeights", "weight_type1 must be -1, 0, 1, 2 or 3 (see options_forecast)."
        )
    try:
        AICc = np.asarray(AICc, dtype=float)
        ok = AICc.size > 0 and (AICc.ndim <= 1 or (AICc.ndim == 2 and 1 in AICc.shape))
        ok = ok and bool(np.all(np.isfinite(AICc)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleWeights", "AICc must hold one finite value per member model."
        )
    AICc = AICc.ravel()
    k = AICc.size
    if weight_type1 == -1:
        weights = np.ones(k) / k
    elif weight_type1 in (0, 1):
        if weight_type1 == 0:
            subepidemicWarning(
                "SubEpiPredict:WeightCodeZero",
                "weight_type1 = 0 used 1/AICc, which is not a weight (equal weights for count "
                "likelihoods, negative weights when AICc changes sign). It now gives Akaike "
                "weights, as weight_type1 = 1 does.",
            )
        relativeLikelihood = np.exp(-(AICc - np.min(AICc)) / 2)
        weights = relativeLikelihood / np.sum(relativeLikelihood)
    elif weight_type1 == 2:
        weights = _inverseScore(WISC, k, "WISC", 2)
    else:
        weights = _inverseScore(WISF, k, "WISF", 3)
    if not (np.all(np.isfinite(weights)) and np.all(weights >= 0) and abs(np.sum(weights) - 1) < 1e-12):
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleWeights",
            "The ensemble weights are not finite, nonnegative numbers that sum to 1.",
        )
    return weights


def _localNumber(v):
    if math.isinf(v):
        return "Inf"
    if v < 1e4:
        return "%.3g" % v
    return "%.2g" % v


def subepidemicEvidenceRatios(AICc, topmodels, area, limit=None):
    """Evidence ratios of the top-ranked models, with a warning for weak ones.

    Port of ``subepidemicEvidenceRatios.m``. ``AICc`` holds the AICc of the
    ranked candidates, best first; ``topmodels`` the number of top-ranked
    models kept; ``area`` the area, for the message; ``limit`` the evidence
    ratio from which a model counts as weakly supported (default 100; inf: no
    warning).

    Returns ``(evidenceRatio, weak)``: ``evidenceRatio[k] = exp((AICc[k] -
    AICc[0])/2)``, how many times more likely ranked model 1 is than ranked
    model k+1; ``weak`` the ONE-BASED ranks with a ratio of at least
    ``limit``. Report only: weak models are not dropped.
    """
    if limit is None:
        limit = 100
    AICc = np.asarray(AICc, dtype=float).ravel()
    topmodels = int(min(max(0, math.floor(float(topmodels))), AICc.size))
    evidenceRatio = np.full(topmodels, np.nan)
    weak = np.zeros(0, dtype=int)
    if topmodels == 0:
        return evidenceRatio, weak
    if not math.isfinite(AICc[0]):
        # no usable best model: there is nothing to compare with
        return evidenceRatio, weak
    delta = AICc[:topmodels] - AICc[0]
    with np.errstate(over="ignore"):
        evidenceRatio = np.exp(delta / 2)  # Inf where AICc is Inf
    evidenceRatio[0] = 1
    weak = np.flatnonzero(evidenceRatio >= limit) + 1
    if weak.size == 0:
        return evidenceRatio, weak
    ratioText = ", ".join(_localNumber(v) for v in evidenceRatio[weak - 1])
    deltaText = ", ".join("%.1f" % v for v in delta[weak - 1])
    if weak.size == 1:
        who = "ranked model %d has an evidence ratio" % weak[0]
        them = "it"
        capital = "It is"
    else:
        who = "ranked models %s have evidence ratios" % ", ".join(str(v) for v in weak)
        them = "them"
        capital = "They are"
    subepidemicWarning(
        "SubEpiPredict:WeakTopModel",
        "Area %d: %s of %s against ranked model 1 (AICc higher by %s). The data give %s little support "
        "(evidence ratio of %s or more). %s kept, because topmodelsx = %d; an ensemble with equal weights "
        "(weight_type1 = -1) counts %s as much as the best model.",
        area, who, ratioText, deltaText, them, _localNumber(limit), capital, topmodels, them,
    )  # fmt: skip
    return evidenceRatio, weak
