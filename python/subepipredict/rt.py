"""Effective reproduction number from an incidence curve (``get_Rt.m``)."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import expon, gamma

from .errors import SubEpiPredictError


def _is_real_scalar(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def _colon(a, step, b):
    """MATLAB's ``a:step:b`` for a positive step (the last point never exceeds b)."""
    n = int(math.floor((b - a) / step + 1e-10)) + 1
    if n <= 0:
        return np.zeros(0)
    return a + step * np.arange(n)


def get_Rt(timevect, incidence1, type_GId1, mean1, var1=None, minCoverage=None):
    """Effective reproduction number from an incidence curve (renewal equation).

    Port of ``get_Rt.m``.

    Parameters
    ----------
    timevect, incidence1 : times and model incidence (same length)
    type_GId1 : generation-interval distribution: 1 gamma, 2 exponential,
        3 delta (fixed interval)
    mean1, var1 : mean and VARIANCE of the generation interval, in the time
        units of ``timevect`` (``var1`` is used by the gamma only)
    minCoverage : default 0.99. Rt is reported only where the curve so far
        covers at least this share of the generation-interval distribution;
        before that it is NaN. 0 turns the blanking off.

    Returns
    -------
    ``(Rs, ts, coverage)``
        ``Rs[k]`` is the estimate at time ``ts[k+1]``; ``ts`` is the grid of
        step 0.01 from the first to the last time; ``coverage[k]`` is the
        share of the generation-interval distribution that lies within the
        curve at ``ts[k+1]``.

    The estimate is ``R(t) = I(t) / sum over past lags s of I(t-s) w(s)``. At a
    time ``t`` after the start of the curve only lags up to ``t`` are
    available, so the denominator misses the share ``1 - coverage`` of the
    infectors and R(t) is too large by a factor of about ``1/coverage``; that
    is why early values are blanked.
    """
    if minCoverage is None or (np.ndim(minCoverage) > 0 and np.size(minCoverage) == 0):
        minCoverage = 0.99
    if not (_is_real_scalar(minCoverage) and 0 <= minCoverage <= 1):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidRtSettings",
            "The minimum generation-interval coverage must be one number from 0 to 1.",
        )
    if not (_is_real_scalar(type_GId1) and type_GId1 in (1, 2, 3)):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidRtSettings", "type_GId1 must be 1 (gamma), 2 (exponential) or 3 (delta)."
        )
    if not (_is_real_scalar(mean1) and math.isfinite(mean1) and mean1 > 0):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidRtSettings",
            "The mean of the generation interval must be one finite positive number.",
        )
    if type_GId1 == 1 and not (_is_real_scalar(var1) and math.isfinite(var1) and var1 > 0):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidRtSettings", "The gamma generation interval needs a finite positive variance."
        )
    timevect = np.asarray(timevect, dtype=float).ravel()
    incidence1 = np.asarray(incidence1, dtype=float).ravel()
    if timevect.size != incidence1.size:
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidRtSettings", "timevect and the incidence curve must have the same length."
        )
    if timevect.size < 2:
        return np.zeros(0), timevect, np.zeros(0)  # one point: nothing to estimate

    dt1 = 0.01
    timevect2 = _colon(timevect[0], dt1, timevect[-1])
    incidence1 = np.interp(timevect2, timevect, incidence1)
    timevect3 = timevect2 - np.min(timevect2)
    ts = timevect2
    n = incidence1.size
    if n < 2:
        return np.zeros(0), ts, np.zeros(0)

    # Generation-interval weights: prob1[j-1] is the probability of a lag in
    # ((j-1)*dt1, j*dt1], used at lag j*dt1.
    if type_GId1 == 1:
        # gamma
        b = var1 / mean1
        a = mean1 / b
        prob1 = np.diff(gamma.cdf(timevect3, a, scale=b))
    elif type_GId1 == 2:
        # exponential
        prob1 = np.diff(expon.cdf(timevect3, scale=mean1))
    else:
        # delta: all weight at lag mean1 (to the nearest grid step, at least one)
        lag1 = max(1, int(round(mean1 / dt1)))
        prob1 = np.zeros(max(n - 1, lag1))
        prob1[lag1 - 1] = 1

    # Renewal denominator by convolution: at time index i it is
    # sum_{j=1}^{i} I(i-j) * prob1(j); only the first n-1 values are needed.
    w = prob1[: n - 1]
    c = np.convolve(incidence1[: n - 1], w, mode="full")
    denominator = c[: n - 1]

    Rs = np.full(n - 1, np.nan)
    currentIncidence = incidence1[1:n]
    valid = denominator > 0
    Rs[valid] = currentIncidence[valid] / denominator[valid]

    # Share of the generation-interval distribution inside the curve at each time.
    coverage = np.zeros(n - 1)
    coverage[: w.size] = np.cumsum(w)
    if w.size < n - 1:
        coverage[w.size :] = coverage[w.size - 1] if w.size else 0.0
    coverage = np.minimum(coverage, 1)

    if minCoverage > 0:
        Rs[coverage < minCoverage - 1e-9] = np.nan
    return Rs, ts, coverage
