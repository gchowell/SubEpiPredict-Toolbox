"""AICc and the number of estimated parameters (``getAICc.m``, ``get_nparams.m``)."""

from __future__ import annotations

import math


def get_nparams(method1, dist1, npatches, flag1, fixI0, onset_fixed):
    """Number of estimated parameters of an n-sub-epidemic model.

    Port of ``get_nparams.m``. Per sub-epidemic: GGM 2, GLM 3, GRM 4,
    logistic 2, Richards 3, Gompertz 2. Plus 1 for the onset threshold of an
    asynchronous model with more than one sub-epidemic, 1 when the initial
    value is estimated (``fixI0 = 0``), and the parameters of the error
    structure: 1 for the normal variance (method 0 with dist 0), 1 for alpha
    (methods 3 and 4), 2 for alpha and d (method 5), 1 for the Laplace scale
    (method 6).
    """
    perPatch = {0: 2, 1: 3, 2: 4, 3: 2, 4: 3, 5: 2}
    numparams = 0
    for _ in range(int(npatches)):
        numparams += perPatch.get(int(flag1), 0)  # flag: 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz

    if npatches > 1 and onset_fixed == 0:
        numparams += 1  # Cthr

    if fixI0 == 0:  # fix initial datum or estimated
        numparams += 1

    if method1 == 0 and dist1 == 0:  # Normal distribution -- one parameter for variance
        numparams += 1
    elif method1 == 3 or method1 == 4:  # Neg. Binomial requires one more parameter (alpha)
        numparams += 1
    elif method1 == 5:
        numparams += 2  # Neg. Binomial requires 2 more parameters (alpha,d)
    elif method1 == 6:  # Laplace distribution for SAD method -- one parameter for the scale parameter
        numparams += 1
    return numparams


def getAICc(method1, dist1, npatches, flag1, fixI0, fval, n, onset_fixed):
    """Corrected Akaike information criterion of a fitted model.

    Port of ``getAICc.m``. ``fval`` is the objective at the optimum (SSE for
    method 0, sum of absolute deviations for method 6, negative log-likelihood
    without its data constant otherwise) and ``n`` the number of data points.

    Returns ``(AICc, part1, part2, numparams)``; ``part1`` is the goodness of
    fit and ``part2`` the penalty. The small-sample term is defined only for
    ``n > numparams + 1``; otherwise ``AICc = inf`` so that a model with more
    parameters than the data support can never be ranked ahead of a valid one.
    """
    numparams = get_nparams(method1, dist1, npatches, flag1, fixI0, onset_fixed)

    if n - numparams - 1 <= 0:
        return math.inf, math.nan, math.inf, numparams

    correction = 2 * numparams + (2 * numparams * (numparams + 1)) / (n - numparams - 1)
    with_log = lambda v: math.log(v) if v > 0 else (-math.inf if v == 0 else math.nan)  # noqa: E731
    if method1 == 0:
        part1 = n * with_log(fval)
        part2 = correction
    elif method1 == 6:
        part1 = 2 * n * with_log(2 * fval / n)
        part2 = 2 * n + correction
    else:
        part1 = -2 * (-fval)
        part2 = correction
    return part1 + part2, part1, part2, numparams
