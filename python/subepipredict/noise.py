"""Error structures: simulated incidence curves and the empirical variance-to-mean ratio.

Ports of ``AddPoissonError.m`` and ``getMeanVarianceRatio.m``.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError, subepidemicWarning
from .utils import get_rng, mat2str


def _show(v):
    try:
        a = np.asarray(v, dtype=float)
    except (TypeError, ValueError):
        return "a " + type(v).__name__
    if a.size == 0:
        return "an empty value"
    if a.size <= 4:
        return mat2str(a.ravel() if a.ndim else a, 6)
    return "a " + type(v).__name__


def _real_scalar(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def AddPoissonError(yi, M, dist1, factor1=None, d=None, rng=None):
    """Simulated incidence curves around a cumulative model curve.

    Port of ``AddPoissonError.m``.

    Parameters
    ----------
    yi : cumulative model curve (a vector)
    M : number of curves to draw
    dist1 : error structure

        0  normal, standard deviation ``factor1``
        1  Poisson
        2  negative binomial, VAR = factor1*mean
        3  negative binomial, VAR = mean + factor1*mean
        4  negative binomial, VAR = mean + factor1*mean^2
        5  negative binomial, VAR = mean + factor1*mean^d
        6  Laplace, scale ``factor1``
    factor1 : noise scale or dispersion (not used for ``dist1 = 1``)
    d : power of the mean for ``dist1 = 5`` (default 1)
    rng : ``numpy.random.Generator``; the toolbox stream when omitted

    Returns
    -------
    ``curves``, ``len(yi)``-by-``M``. Its first row is ``yi[0]`` in every
    column: the first observation is not perturbed. Row ``t+1`` is drawn
    around the increment ``yi[t+1] - yi[t]``. A zero expected count gives 0;
    zero overdispersion gives the Poisson limit; where the cumulative curve
    decreases the expected count is taken as 0. The result never contains NaN
    or Inf.
    """
    rng = get_rng(rng)
    if d is None or (np.ndim(d) > 0 and np.size(d) == 0):
        d = 1
    try:
        yi = np.asarray(yi, dtype=float)
        ok = yi.size > 0 and (yi.ndim <= 1 or (yi.ndim == 2 and 1 in yi.shape)) and bool(np.all(np.isfinite(yi)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:ErrorGenerator", "The cumulative curve must be a non-empty vector of finite real numbers."
        )
    yi = yi.ravel()
    if not (_real_scalar(M) and math.isfinite(M) and M >= 0 and M == int(M)):
        raise SubEpiPredictError("SubEpiPredict:ErrorGenerator", "The number of curves must be a nonnegative integer.")
    M = int(M)
    if not (_real_scalar(dist1) and dist1 in (0, 1, 2, 3, 4, 5, 6)):
        raise SubEpiPredictError("SubEpiPredict:ErrorGenerator", "dist1 must be an integer from 0 to 6.")
    if dist1 != 1:
        if not (_real_scalar(factor1) and math.isfinite(factor1) and factor1 >= 0):
            raise SubEpiPredictError(
                "SubEpiPredict:ErrorGenerator",
                "dist1 = %d needs factor1, the noise scale or dispersion, as one finite nonnegative "
                "number; received %s. Check that method1 and dist1 in the options match." % (dist1, _show(factor1)),
            )
        factor1 = float(factor1)
    if dist1 == 5 and not (_real_scalar(d) and math.isfinite(d) and d > 0):
        raise SubEpiPredictError("SubEpiPredict:ErrorGenerator", "dist1 = 5 needs a finite positive power d.")
    if dist1 == 2 and factor1 < 1:
        subepidemicWarning(
            "SubEpiPredict:UnderdispersedFactor",
            "dist1 = 2 with a variance-to-mean factor of %g (below 1) is not a negative binomial. "
            "Poisson noise (factor 1) is used instead.",
            factor1,
        )
        factor1 = 1.0

    mu = np.diff(yi)
    if dist1 in (1, 2, 3, 4, 5):
        # Count means cannot be negative: a decreasing stretch has expected count 0.
        # ODE noise on a flat curve is passed over; a real decrease is reported.
        tolerance = 1e-6 * max(1.0, float(np.max(np.abs(yi))))
        if np.any(mu < -tolerance):
            subepidemicWarning(
                "SubEpiPredict:DecreasingCurve",
                "The cumulative model curve decreases at %d of %d steps (largest drop %g). "
                "The expected count is taken as 0 there.",
                int(np.sum(mu < -tolerance)), mu.size, float(np.max(-mu)),
            )  # fmt: skip
        mu = np.maximum(mu, 0)

    curves = np.zeros((yi.size, M))
    curves[0, :] = yi[0]
    if M == 0 or mu.size == 0:
        return curves

    if dist1 == 0:
        # Normal distribution
        curves[1:, :] = mu[:, None] + factor1 * rng.standard_normal((mu.size, M))
    elif dist1 == 6:
        # Laplace distribution around abs(increment), by its inverse CDF
        eps = np.finfo(float).eps
        U = np.minimum(np.maximum(rng.random((mu.size, M)), eps), 1 - eps) - 0.5
        curves[1:, :] = np.abs(mu)[:, None] - factor1 * np.sign(U) * np.log1p(-2 * np.abs(U))
    else:
        # Poisson and negative binomial: VAR = mean*(1 + extraRatio)
        for t in range(mu.size):
            mean1 = float(mu[t])
            if mean1 == 0:
                continue  # a zero expected count is exactly 0
            try:
                if dist1 == 1:
                    extraRatio = 0.0
                elif dist1 == 2:
                    extraRatio = factor1 - 1  # VAR = factor1*mean
                elif dist1 == 3:
                    extraRatio = factor1  # VAR = mean + alpha*mean
                elif dist1 == 4:
                    extraRatio = factor1 * mean1  # VAR = mean + alpha*mean^2
                else:
                    extraRatio = factor1 * mean1 ** (d - 1)  # VAR = mean + alpha*mean^d
            except OverflowError:
                extraRatio = math.inf
            if not math.isfinite(extraRatio) or extraRatio < 0:
                raise SubEpiPredictError(
                    "SubEpiPredict:ErrorGenerator", "The variance of the error structure is invalid or overflowed."
                )
            p1 = 1 / (1 + extraRatio)
            if extraRatio == 0 or p1 == 1:
                curves[t + 1, :] = rng.poisson(mean1, M)
            else:
                r1 = mean1 / extraRatio
                if not math.isfinite(r1) or r1 <= 0 or p1 <= 0:
                    raise SubEpiPredictError(
                        "SubEpiPredict:ErrorGenerator",
                        "The negative binomial parameters are outside their numerical range.",
                    )
                # Gamma-Poisson mixture: the same distribution as nbinrnd(r1, p1),
                # valid for any positive r1.
                lam = rng.gamma(r1, (1 - p1) / p1, M)
                curves[t + 1, :] = rng.poisson(lam)

    if not np.all(np.isfinite(curves)):
        raise SubEpiPredictError(
            "SubEpiPredict:ErrorGenerator", "The error generator produced a value that is not a finite real number."
        )
    return curves


def getMeanVarianceRatio(data1, binsize1, flag1):
    """Variance-to-mean ratio of the data in windows of ``binsize1`` points.

    Port of ``getMeanVarianceRatio.m``. ``flag1 = 1``: consecutive windows;
    ``flag1 = 2``: a window starting at every point. Returns
    ``(ratios, timebins1)``; the columns of ``ratios`` are the ratio, the mean
    and the (sample) variance of each window.
    """
    data1 = np.asarray(data1, dtype=float).ravel()
    binsize1 = int(binsize1)
    ratios = []
    timebins1 = []
    if flag1 == 1:  # no moving window
        starts = range(0, data1.size - binsize1, binsize1)
    elif flag1 == 2:  # moving window
        starts = range(0, data1.size - binsize1, 1)
    else:
        starts = []
    with np.errstate(all="ignore"):
        for i in starts:
            window = data1[i : i + binsize1]
            mean1 = float(np.mean(window))
            var1 = float(np.var(window, ddof=1)) if window.size > 1 else 0.0
            ratios.append([np.divide(var1, mean1), mean1, var1])
            timebins1.append(i + (binsize1 + 1) / 2 if flag1 == 1 else i + 1)
    return np.array(ratios, dtype=float).reshape(-1, 3), np.array(timebins1, dtype=float)
