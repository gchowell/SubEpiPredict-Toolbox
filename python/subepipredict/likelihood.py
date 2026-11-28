"""Fitting objective and the negative-binomial likelihood.

Ports of ``subepidemicNegativeBinomialNLL.m`` and ``subepidemicObjective.m``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln

from .errors import SubEpiPredictError
from .model import simulateSubepidemic
from .utils import Struct

_LOG16 = math.log(16.0)


def _is_real_scalar(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def subepidemicNegativeBinomialNLL(y, mu, alpha, d, method):
    """Negative-binomial fitting objective of the toolbox.

    Port of ``subepidemicNegativeBinomialNLL.m``.

    ``method`` 3: variance = mu + alpha*mu; 4: variance = mu + alpha*mu^2;
    5: variance = mu + alpha*mu^d.

    The value is the negative log-likelihood WITHOUT ``sum(gammaln(y+1))``, as
    in the original objective. For fractional ``y`` the legacy convention is
    kept: the product term uses ``floor(y)`` factors,
    ``gammaln(r + floor(y)) - gammaln(r)``, and the other terms ``y`` itself.
    ``mu`` must be strictly positive. Invalid numeric values return ``inf``;
    shape and method errors raise. At ``alpha = 0`` integer ``y`` give the
    Poisson limit and fractional ``y`` give ``inf``.
    """
    try:
        y = np.asarray(y)
        mu = np.asarray(mu)
        shape_ok = (
            y.dtype.kind in "biufc"
            and mu.dtype.kind in "biufc"
            and y.size > 0
            and mu.size == y.size
            and (y.ndim <= 1 or (y.ndim == 2 and 1 in y.shape))
            and (mu.ndim <= 1 or (mu.ndim == 2 and 1 in mu.shape))
        )
    except (TypeError, ValueError):
        shape_ok = False
    if not shape_ok:
        raise SubEpiPredictError(
            "SubEpiPredict:NBInputShape", "Y and MU must be nonempty numeric vectors of equal length."
        )
    if not (_is_real_scalar(method) and method in (3, 4, 5)):
        raise SubEpiPredictError("SubEpiPredict:NBMethod", "NB method must be 3, 4, or 5.")

    nll = math.inf
    if (
        y.dtype.kind == "c"
        or mu.dtype.kind == "c"
        or not (_is_real_scalar(alpha) and math.isfinite(alpha) and alpha >= 0)
        or not (_is_real_scalar(d) and math.isfinite(d))
    ):
        return nll
    y = y.astype(float).ravel()
    mu = mu.astype(float).ravel()
    if np.any(~np.isfinite(y) | (y < 0) | ~np.isfinite(mu) | (mu <= 0)):
        return nll
    alpha = float(alpha)
    d = float(d)
    n = np.floor(y)
    logMu = np.log(mu)

    if alpha == 0:
        if np.all(y == n):
            nll = float(np.sum(mu - y * logMu))
            if not math.isfinite(nll):
                nll = math.inf
        return nll

    logAlpha = math.log(alpha)
    if method == 3:
        logBeta = np.full(mu.shape, logAlpha)
        logR = logMu - logAlpha
    elif method == 4:
        logBeta = logAlpha + logMu
        logR = np.full(mu.shape, -logAlpha)
    else:
        logBeta = logAlpha + (d - 1) * logMu
        logR = -logAlpha + (2 - d) * logMu
    if np.any(~np.isfinite(logBeta) | ~np.isfinite(logR)):
        return nll

    with np.errstate(all="ignore"):
        # beta = mu/r. Stable log(1+beta), log(beta/(1+beta)), and r*log(1+beta).
        logOnePlusBeta = np.maximum(logBeta, 0) + np.log1p(np.exp(-np.abs(logBeta)))
        logFailure = np.zeros(y.shape)
        penalty = np.zeros(y.shape)
        smallBeta = logBeta < 0
        if np.any(smallBeta):
            b = np.exp(logBeta[smallBeta])
            ratio = np.ones(b.shape)
            nonzero = b > 0
            ratio[nonzero] = np.log1p(b[nonzero]) / b[nonzero]
            logFailure[smallBeta] = logBeta[smallBeta] - np.log1p(b)
            penalty[smallBeta] = mu[smallBeta] * ratio
        big = ~smallBeta
        if np.any(big):
            logFailure[big] = -np.log1p(np.exp(-logBeta[big]))
            penalty[big] = np.exp(logR[big] + np.log(logOnePlusBeta[big]))

        logTerms = np.zeros(y.shape)
        largeR = logR >= _LOG16
        if np.any(~largeR):
            # Gamma recurrence removes gammaln(r)'s singularity as r approaches 0:
            # log Gamma(r+n)-log Gamma(r) = log(r)+log Gamma(r+n)-log Gamma(r+1).
            ii = np.flatnonzero(~largeR)
            ni = n[ii]
            r = np.exp(logR[ii])
            rising = np.zeros(ni.shape)
            positive = ni > 0
            rising[positive] = logR[ii[positive]] + gammaln(r[positive] + ni[positive]) - gammaln(1 + r[positive])
            logTerms[ii] = rising + y[ii] * logFailure[ii] - penalty[ii]
        if np.any(largeR):
            ii = np.flatnonzero(largeR)
            ni = n[ii]
            correction = _scaledLogRisingFactorial(ni, logR[ii])
            # Cancel n*log(r) against n*log(beta) analytically, not numerically.
            # The (y-n)*log(beta) term MUST remain for legacy fractional inputs.
            logTerms[ii] = (
                ni * logMu[ii] + (y[ii] - ni) * logBeta[ii] + correction - y[ii] * logOnePlusBeta[ii] - penalty[ii]
            )
    if np.all(np.isfinite(logTerms)):
        nll = float(-np.sum(logTerms))
        if not math.isfinite(nll):
            nll = math.inf
    return nll


def _scaledLogRisingFactorial(n, logR):
    """log Gamma(r+n) - log Gamma(r) - n*log(r), for r >= 16 and integer n >= 0.

    With t = n/r and H(t) = ((1+t)*log1p(t) - t)/t the Stirling main term is
    n*H(t) - 0.5*log1p(t), which avoids subtracting quantities of order
    r*log(r). Corrections through r^(-9).
    """
    correction = np.zeros(n.shape)
    active = n > 1  # n=0 and n=1 give exactly zero, including at r=Inf.
    if not np.any(active):
        return correction
    nn = n[active]
    u = np.exp(-logR[active])
    t = nn * u
    L = np.log1p(t)
    H = np.zeros(t.shape)
    small = t <= 0.01
    s = t[small]
    # Taylor polynomial through t^9; first omitted term <= t^10/110.
    H[small] = s * (
        1 / 2
        + s * (-1 / 6 + s * (1 / 12 + s * (-1 / 20 + s * (1 / 30 + s * (-1 / 42 + s * (1 / 56 + s * (-1 / 72 + s / 90)))))))
    )
    H[~small] = (1 + 1 / t[~small]) * L[~small] - 1
    # S(r+n)-S(r); expm1 retains accuracy when n/r is tiny.
    u2 = u * u
    stirlingDifference = u * (
        np.expm1(-L) / 12
        + u2 * (-np.expm1(-3 * L) / 360 + u2 * (np.expm1(-5 * L) / 1260 + u2 * (-np.expm1(-7 * L) / 1680 + u2 * np.expm1(-9 * L) / 1188)))
    )
    correction[active] = nn * H - 0.5 * L + stirlingDifference
    return correction


_CONTEXT_FIELDS = (
    "timevect", "ydata", "I0", "npatches", "onset_fixed", "onset_thr", "flag1", "method1", "odeOptions",
)  # fmt: skip


def subepidemicObjective(z, context):
    """Fitting objective with an explicit context.

    Port of ``subepidemicObjective.m``.

    Parameters
    ----------
    z : parameter vector ``[r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]``
    context : dict / Struct with the fields ``timevect``, ``ydata``, ``I0``,
        ``npatches``, ``onset_fixed``, ``onset_thr``, ``flag1``, ``method1``
        (0 least squares, 1 Poisson, 3/4/5 negative binomial, 6 sum of
        absolute deviations) and ``odeOptions``.

    Returns
    -------
    ``(value, detail)``
        ``value`` is the SSE, the negative log-likelihood without its
        data-only constant, or the sum of absolute deviations. ``detail`` is a
        Struct: ``success``, ``message``, ``nActive`` (sub-epidemics active by
        the end of the grid), ``activation``, ``curve`` (fitted incidence).

    A parameter vector is invalid (``detail.success`` False, ``value`` 1e10)
    when the simulation fails, the trajectory is identically zero, or the
    objective is not a finite real number.
    """
    value = 1e10
    detail = Struct(success=False, message="", nActive=0, activation=None, curve=None)
    try:
        n = int(context["npatches"])
        try:
            zz = np.asarray(z, dtype=float).ravel()
        except (TypeError, ValueError):
            zz = np.empty(0)
        if not (zz.size == 4 * n + 2 and np.all(np.isfinite(zz))):
            detail.message = "Parameter vector has the wrong size or nonfinite values."
            return value, detail
        _, _, yfit, _, activation = simulateSubepidemic(
            context["timevect"], context["I0"], n, context["onset_fixed"], context["onset_thr"], context["flag1"],
            zz[0:n], zz[n : 2 * n], zz[2 * n : 3 * n], zz[3 * n : 4 * n], context["odeOptions"],
        )  # fmt: skip
        alpha = float(zz[-2])
        d = float(zz[-1])
        detail.activation = activation
        detail.nActive = int(np.sum(activation.active))
        if np.sum(yfit) == 0:
            detail.message = "The fitted trajectory is identically zero."
            return value, detail
        yfit = np.where(yfit == 0, 0.001, yfit)  # legacy convention, kept so objectives are comparable
        ydata = np.asarray(context["ydata"], dtype=float).ravel()
        method1 = context["method1"]
        with np.errstate(all="ignore"):
            if method1 == 0:
                v = float(np.sum((ydata - yfit) ** 2))
            elif method1 == 1:
                v = float(-np.sum(ydata * np.log(yfit) - yfit))
            elif method1 in (3, 4, 5):
                v = subepidemicNegativeBinomialNLL(ydata, yfit, alpha, d, method1)
            elif method1 == 6:
                v = float(np.sum(np.abs(ydata - yfit)))
            else:
                raise SubEpiPredictError("SubEpiPredict:InvalidMethod", "method1 must be 0, 1, 3, 4, 5 or 6.")
        if not math.isfinite(v):
            detail.message = "The objective is not a finite real number at this parameter vector."
            return value, detail
        value = v
        detail.success = True
        detail.curve = yfit
    except SubEpiPredictError as ME:
        value = 1e10
        detail.success = False
        detail.message = "%s: %s" % (ME.identifier, ME.message)
    except Exception as ME:  # any other failure also marks the vector as invalid
        value = 1e10
        detail.success = False
        detail.message = "%s: %s" % (type(ME).__name__, ME)
    return value, detail
