"""Parameter identifiability: PII, bound hits, Monte Carlo error, reports.

Ports of ``subepidemicIdentifiability.m``, ``subepidemicPII.m``,
``subepidemicBoundHitFractions.m``, ``subepidemicBoundSummary.m``,
``subepidemicQuantileMCSE.m``, ``sci_interval.m``, ``subepidemicPrintPII.m``
and ``subepidemicWritePII.m``.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError
from .search import subepidemicUnitMap
from .utils import Struct, num2str, quantile


def sci_interval(lb, ub):
    """log10 width of a 95% interval, the toolbox's identifiability score.

    Port of ``sci_interval.m``: ``log10(UB/LB)`` with the lower limit floored
    at ``1e-6*UB``. ``LB == UB`` gives 0; ``LB == 0 < UB`` gives 6, the cap of
    the score ("not bounded away from zero"); NaN in either limit gives NaN.
    """
    lb = np.asarray(lb, dtype=float)
    ub = np.asarray(ub, dtype=float)
    floorRatio = 1e-6
    with np.errstate(all="ignore"):
        s = np.log10(ub / np.fmax(lb, floorRatio * ub))
    s = np.where(ub == lb, 0.0, s)
    s = np.where(np.isnan(lb) | np.isnan(ub), np.nan, s)
    return s if s.ndim else float(s)


def subepidemicPII(lo, hi, LBbox, UBbox, scaleUser=None):
    """Parameter identifiability index (PII) and its companions.

    Port of ``subepidemicPII.m``. ``lo``, ``hi`` are the 2.5% and 97.5%
    bootstrap quantiles; ``LBbox``, ``UBbox`` the search bounds.

    Returns ``(PII, rangeFrac, scale, PIIbox, ratio)``:

    ``PII``        ``log10(hi/lo)``, the width of the 95% interval in decades;
                   NaN unless ``lo > 0``
    ``rangeFrac``  ``(hi - lo)/(UBbox - LBbox)``, in [0, 1]; NaN for a fixed
                   parameter
    ``scale``      headline measure: 1 log scale (lower bound positive), 0
                   range fraction; ``scaleUser`` overrides it per parameter
                   (1, 0, or NaN for automatic)
    ``PIIbox``     ``log10(UBbox/LBbox)``, the largest PII the bounds allow
    ``ratio``      ``PII/PIIbox`` on the log scale, ``rangeFrac`` otherwise;
                   near 1 = the interval is the search range
    """
    lo = np.asarray(lo, dtype=float).ravel()
    hi = np.asarray(hi, dtype=float).ravel()
    LBbox = np.asarray(LBbox, dtype=float).ravel()
    UBbox = np.asarray(UBbox, dtype=float).ravel()
    n = lo.size
    if hi.size != n or LBbox.size != n or UBbox.size != n:
        raise SubEpiPredictError(
            "SubEpiPredict:PIISize", "lo, hi, LBbox and UBbox must have the same number of elements."
        )

    with np.errstate(all="ignore"):
        # Log-scale index: only where both quantiles are positive.
        PII = np.full(n, np.nan)
        ok = (lo > 0) & (hi > 0)
        PII[ok] = np.log10(hi[ok]) - np.log10(lo[ok])

        # Share of the search range covered by the interval.
        width = UBbox - LBbox
        rangeFrac = np.full(n, np.nan)
        free = (width > 0) & np.isfinite(width) & np.isfinite(lo) & np.isfinite(hi)
        rangeFrac[free] = np.minimum(np.maximum((hi[free] - lo[free]) / width[free], 0), 1)

        # Headline scale: automatic rule, then the user's override where given.
        scale = (LBbox > 0).astype(float)
        if scaleUser is not None and np.size(scaleUser) > 0:
            scaleUser = np.asarray(scaleUser, dtype=float).ravel()
            m = min(scaleUser.size, n)
            padded = np.full(n, np.nan)
            padded[:m] = scaleUser[:m]
            given = np.isfinite(padded)
            if np.any((padded[given] != 0) & (padded[given] != 1)):
                raise SubEpiPredictError(
                    "SubEpiPredict:PIIScale",
                    "PIIscale entries must be 1 (log scale), 0 (range fraction) or NaN (automatic).",
                )
            scale[given] = padded[given]

        # Largest PII the bounds allow; undefined for a nonpositive lower bound.
        PIIbox = np.full(n, np.nan)
        pos = (LBbox > 0) & (UBbox > 0)
        PIIbox[pos] = np.log10(UBbox[pos]) - np.log10(LBbox[pos])

        # Bound-censoring diagnostic on the headline scale of each parameter.
        ratio = np.full(n, np.nan)
        isLog = scale == 1
        ratio[isLog] = PII[isLog] / PIIbox[isLog]
        ratio[~isLog] = rangeFrac[~isLog]
        ratio[isLog & (PIIbox == 0)] = np.nan  # fixed parameter: 0/0, made explicit
    return PII, rangeFrac, scale, PIIbox, ratio


def subepidemicBoundHitFractions(Phatss, LB, UB, relTol=None):
    """Share of bootstrap refits sitting on the search bounds.

    Port of ``subepidemicBoundHitFractions.m``. A refit counts as hitting the
    lower bound when ``P <= LB + relTol*|LB|`` and the upper bound when
    ``P >= UB - relTol*|UB|`` (``relTol`` default 1e-3; for a bound that is
    exactly 0 the hit must be exact).

    Returns ``(fracLB, fracUB, hitLB, hitUB)``: the fractions per column (NaN
    for fixed columns and for columns whose refits are all NaN) and the
    M-by-d masks of the individual hits.
    """
    if relTol is None:
        relTol = 1e-3
    Phatss = np.asarray(Phatss, dtype=float)
    if Phatss.ndim == 1:
        Phatss = Phatss.reshape(1, -1)
    M, d = Phatss.shape
    LB = np.asarray(LB, dtype=float)
    UB = np.asarray(UB, dtype=float)

    def broadcast(B):
        if B.ndim <= 1 or 1 in B.shape:
            if not (d == 1 and B.ndim == 2 and B.shape[0] == M and M > 1):
                B = np.tile(B.ravel(), (M, 1))
        return B

    LB = broadcast(LB)
    UB = broadcast(UB)
    if LB.shape[0] != M or UB.shape[0] != M or LB.shape[1] < d or UB.shape[1] < d:
        raise SubEpiPredictError(
            "SubEpiPredict:BoundHitSizeMismatch",
            "LB/UB must be 1-by-d (or longer) or M-by-d (or wider) to match Phatss (%d-by-%d)." % (M, d),
        )
    LB = LB[:, :d].copy()
    UB = UB[:, :d].copy()
    # Guard against flipped boxes.
    flip = LB > UB
    if np.any(flip):
        LB[flip], UB[flip] = UB[flip], LB[flip]

    tolLB = relTol * np.abs(LB)
    tolUB = relTol * np.abs(UB)
    valid = ~np.isnan(Phatss)
    with np.errstate(invalid="ignore"):
        hitLB = (Phatss <= LB + tolLB) & valid
        hitUB = (Phatss >= UB - tolUB) & valid
    nValid = valid.sum(axis=0)
    fracLB = hitLB.sum(axis=0) / np.maximum(nValid, 1)
    fracUB = hitUB.sum(axis=0) / np.maximum(nValid, 1)
    # Fixed columns (LB==UB in every refit) and all-NaN columns carry no information.
    pinned = np.all(LB == UB, axis=0)
    fracLB[pinned | (nValid == 0)] = np.nan
    fracUB[pinned | (nValid == 0)] = np.nan
    return fracLB, fracUB, hitLB, hitUB


def subepidemicBoundSummary(P, LB, UB, isLog, npatches):
    """Parameters that sit at a bound, and in how many fits.

    Port of ``subepidemicBoundSummary.m``. ``P`` holds one fitted parameter
    vector per row. A value within 0.1% of a bound counts as at the bound (the
    rule of :func:`subepidemicBoundHitFractions`). Fixed parameters are never
    reported; rows with a non-finite entry are ignored.

    Returns ``(report, lines, nRowsAtBound)``. ``report`` is a list of
    Structs, one per parameter and side with at least one fit at the bound:
    ``name`` ('r1', 'p2', 'K1', 'alpha', 'd'), ``index`` (ZERO-BASED column
    of P), ``side`` ('lower' or 'upper'), ``limit``, ``count``, ``total``,
    ``fraction``, ``setting`` (the fit setting that moves the bound, '' for
    alpha and d) and ``isSearchLimit`` (True for r and K, whose bounds only
    limit the search). ``lines`` holds one printable line per element.
    """
    ok = isinstance(npatches, (int, float, np.integer, np.floating)) and not isinstance(npatches, bool)
    if not (ok and math.isfinite(npatches) and npatches >= 1 and npatches == int(npatches)):
        raise SubEpiPredictError("SubEpiPredict:InvalidPatchCount", "npatches must be a positive integer.")
    npatches = int(npatches)
    nz = 4 * npatches + 2
    LB = np.asarray(LB, dtype=float).ravel()
    UB = np.asarray(UB, dtype=float).ravel()
    isLog = np.asarray(isLog).astype(bool).ravel()
    if LB.size != nz or UB.size != nz or isLog.size != nz:
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Bounds and scale flags must have %d entries for %d sub-epidemics." % (nz, npatches),
        )
    P = np.asarray(P, dtype=float)
    if P.ndim == 1:
        P = P.reshape(1, -1)
    if P.shape[1] != nz:
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Each row must hold %d parameters; received %d columns." % (nz, P.shape[1]),
        )
    P = P[np.all(np.isfinite(P), axis=1), :]
    total = P.shape[0]
    report = []
    lines = []
    if total == 0:
        return report, lines, 0
    m = subepidemicUnitMap("make", dict(z0=LB, LB=LB, UB=UB, isLog=isLog))
    _, _, hitLB, hitUB = subepidemicBoundHitFractions(P, LB, UB)
    nRowsAtBound = int(np.sum(np.any(hitLB[:, m.free] | hitUB[:, m.free], axis=1)))
    letters = "rpaK"
    for j in m.free:
        j = int(j)
        if j < 4 * npatches:
            block = j // npatches
            name = "%s%d" % (letters[block], j - block * npatches + 1)
            prefix = letters[block]
            searchLimit = block == 0 or block == 3
        elif j == 4 * npatches:
            name, prefix, searchLimit = "alpha", "", False
        else:
            name, prefix, searchLimit = "d", "", False
        counts = (int(hitLB[:, j].sum()), int(hitUB[:, j].sum()))
        for s, (side, limit, suffix) in enumerate((("lower", LB[j], "LB"), ("upper", UB[j], "UB"))):
            if counts[s] > 0:
                report.append(
                    Struct(
                        name=name, index=j, side=side, limit=float(limit), count=counts[s], total=total,
                        fraction=counts[s] / total, setting=(prefix + suffix) if prefix else "",
                        isSearchLimit=searchLimit,
                    )  # fmt: skip
                )
                lines.append(
                    "%s at its %s limit (%s) in %d of %d" % (name, side, num2str(limit, 6), counts[s], total)
                )
    return report, lines, nRowsAtBound


def subepidemicQuantileMCSE(Phatss, B2=None):
    """Monte Carlo standard error of the bootstrap interval and PII.

    Port of ``subepidemicQuantileMCSE.m``. The rows of ``Phatss`` are
    resampled with replacement ``B2`` times (default 500); the 2.5% and 97.5%
    quantiles and ``PII = log10(q97.5/q2.5)`` of every column are recomputed
    and the SD across resamples is returned. The resampling uses its own
    random stream with a fixed seed (20261030), so the result is repeatable
    and the toolbox stream is not advanced.

    Returns ``(mcseLB, mcseUB, mcsePII)``. Columns with zero spread and
    all-NaN columns give NaN.
    """
    if B2 is None:
        B2 = 500
    B2 = int(B2)
    Phatss = np.asarray(Phatss, dtype=float)
    if Phatss.ndim == 1:
        Phatss = Phatss.reshape(1, -1)
    M, d = Phatss.shape
    mcseLB = np.full(d, np.nan)
    mcseUB = np.full(d, np.nan)
    mcsePII = np.full(d, np.nan)
    if M < 2 or B2 < 2:
        return mcseLB, mcseUB, mcsePII
    rs = np.random.Generator(np.random.MT19937(20261030))  # local stream: the toolbox stream is untouched
    idx = rs.integers(0, M, size=(M, B2))  # one resample per column
    for k in range(d):
        x = Phatss[:, k]
        if np.all(np.isnan(x)) or np.nanmax(x) - np.nanmin(x) == 0:
            continue  # fixed / degenerate
        xr = x[idx]  # M-by-B2 resampled estimates
        q = quantile(xr, [0.025, 0.975], axis=0)  # 2-by-B2
        mcseLB[k] = np.std(q[0, :], ddof=1)
        mcseUB[k] = np.std(q[1, :], ddof=1)
        ok = (q[0, :] > 0) & (q[1, :] > 0)
        if np.count_nonzero(ok) >= 2:
            mcsePII[k] = np.std(np.log10(q[1, ok]) - np.log10(q[0, ok]), ddof=1)
    return mcseLB, mcseUB, mcsePII


def subepidemicIdentifiability(Phatss, LB, UB, npatches, P0=None, scaleByType=None, B2=None, boundHitTol=None):
    """PII metrics of one ranked sub-epidemic model.

    Port of ``subepidemicIdentifiability.m``.

    ``Phatss``       M-by-(4n+2) bootstrap refits
    ``LB``, ``UB``   search bounds in the same layout
    ``npatches``     n
    ``P0``           the ranked fit to the observed data (None = NaN)
    ``scaleByType``  optional headline scale per parameter type
                     ``[r, p, a, K, alpha, d]``: 1 log scale, 0 range fraction,
                     NaN automatic (setting ``PIIscale``)
    ``B2``           resamples of the Monte Carlo error (default 500)
    ``boundHitTol``  relative tolerance of the bound-hit fractions (1e-3)

    Returns a Struct with one entry per reported parameter (r, p, a, K of
    every sub-epidemic; alpha and d only when estimated): ``name``, ``type``,
    ``subepidemic``, ``column`` (ZERO-BASED column of Phatss), ``estimate``,
    ``median``, ``CI_LB``, ``CI_UB``, ``PII``, ``PII_rangeFrac``,
    ``PII_scale``, ``PII_box``, ``PII_ratio``, ``fracLB``, ``fracUB``,
    ``MCSE_CI_LB``, ``MCSE_CI_UB``, ``MCSE_PII``, ``bound_LB``, ``bound_UB``,
    ``fixed``; and the scalars ``M``, ``B2`` and ``boundHitTol``.

    Reading the numbers: PII near 0 is a tight interval, PII above 1 an
    interval wider than a factor of 10. An interval is evidence about the data
    only when ``PII_ratio`` is well below 1 and ``fracLB``, ``fracUB`` are
    near 0; otherwise the search bounds set its limits.
    """
    ok = isinstance(npatches, (int, float, np.integer, np.floating)) and not isinstance(npatches, bool)
    if not (ok and math.isfinite(npatches) and npatches >= 1 and npatches == int(npatches)):
        raise SubEpiPredictError("SubEpiPredict:InvalidPatchCount", "npatches must be a positive integer.")
    n = int(npatches)
    nz = 4 * n + 2
    Phatss = np.asarray(Phatss, dtype=float)
    if Phatss.ndim == 1:
        Phatss = Phatss.reshape(1, -1)
    LB = np.asarray(LB, dtype=float).ravel()
    UB = np.asarray(UB, dtype=float).ravel()
    if Phatss.shape[1] != nz or LB.size != nz or UB.size != nz:
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Phatss, LB and UB must have %d columns for %d sub-epidemics." % (nz, n),
        )
    if P0 is None or np.size(P0) == 0:
        P0 = np.full(nz, np.nan)
    P0 = np.asarray(P0, dtype=float).ravel()
    if P0.size != nz:
        raise SubEpiPredictError("SubEpiPredict:ParameterDimensionMismatch", "P0 must have %d entries." % nz)
    if scaleByType is None or np.size(scaleByType) == 0:
        scaleByType = np.full(6, np.nan)
    scaleByType = np.asarray(scaleByType, dtype=float).ravel()
    if scaleByType.size != 6:
        raise SubEpiPredictError(
            "SubEpiPredict:PIIScale", "PIIscale must have 6 entries, one per parameter type: [r p a K alpha d]."
        )
    if B2 is None:
        B2 = 500
    if boundHitTol is None:
        boundHitTol = 1e-3

    # One scale entry per column of Phatss.
    scaleUser = np.concatenate([np.repeat(scaleByType[:4], n), scaleByType[4:6]])

    lo = quantile(Phatss, 0.025, axis=0)
    hi = quantile(Phatss, 0.975, axis=0)
    med = np.median(Phatss, axis=0)
    PII, rangeFrac, scale, PIIbox, ratio = subepidemicPII(lo, hi, LB, UB, scaleUser)
    fracLB, fracUB, _, _ = subepidemicBoundHitFractions(Phatss, LB, UB, boundHitTol)
    mcseLB, mcseUB, mcsePII = subepidemicQuantileMCSE(Phatss, B2)

    # Parameters reported: r, p, a, K of every sub-epidemic; alpha and d when estimated.
    keep = list(range(4 * n)) + [4 * n + i for i in range(2) if UB[4 * n + i] > LB[4 * n + i]]
    keep = np.array(keep, dtype=int)
    letters = "rpaK"
    names, types = [], []
    sub = np.full(nz, np.nan)
    for j in range(4 * n):
        block = j // n
        sub[j] = j - block * n + 1
        types.append(letters[block])
        names.append("%s_%d" % (letters[block], sub[j]))
    names += ["alpha", "d"]
    types += ["alpha", "d"]

    pii = Struct()
    pii.name = [names[k] for k in keep]
    pii.type = [types[k] for k in keep]
    pii.subepidemic = sub[keep]
    pii.column = keep
    pii.estimate = P0[keep]
    pii.median = np.asarray(med)[keep]
    pii.CI_LB = np.asarray(lo)[keep]
    pii.CI_UB = np.asarray(hi)[keep]
    pii.PII = PII[keep]
    pii.PII_rangeFrac = rangeFrac[keep]
    pii.PII_scale = scale[keep]
    pii.PII_box = PIIbox[keep]
    pii.PII_ratio = ratio[keep]
    pii.fracLB = fracLB[keep]
    pii.fracUB = fracUB[keep]
    pii.MCSE_CI_LB = mcseLB[keep]
    pii.MCSE_CI_UB = mcseUB[keep]
    pii.MCSE_PII = mcsePII[keep]
    pii.bound_LB = LB[keep]
    pii.bound_UB = UB[keep]
    pii.fixed = LB[keep] == UB[keep]
    pii.M = Phatss.shape[0]
    pii.B2 = B2
    pii.boundHitTol = boundHitTol
    return pii


def subepidemicPrintPII(pii, heading=None, file=None):
    """Print the PII metrics of one ranked model as a table (``subepidemicPrintPII.m``)."""
    if not heading:
        heading = "Parameter identifiability (PII)"
    print("\n%s (%d bootstrap refits):" % (heading, pii.M), file=file)
    print(
        "  %-7s %10s %10s %10s %7s %7s %6s %7s %7s %6s %6s %7s"
        % ("param", "estimate", "CI LB", "CI UB", "PII", "range", "scale", "box", "ratio", "atLB", "atUB", "MCSE"),
        file=file,
    )
    for k in range(len(pii.name)):
        if pii.fixed[k]:
            print("  %-7s %10.4g %s" % (pii.name[k], pii.estimate[k], "  fixed"), file=file)
            continue
        scaleName = "log" if pii.PII_scale[k] == 1 else "range"
        print(
            "  %-7s %10.4g %10.4g %10.4g %7.3f %7.3f %6s %7.3f %7.3f %6.3f %6.3f %7.3f"
            % (pii.name[k], pii.estimate[k], pii.CI_LB[k], pii.CI_UB[k], pii.PII[k], pii.PII_rangeFrac[k], scaleName,
               pii.PII_box[k], pii.PII_ratio[k], pii.fracLB[k], pii.fracUB[k], pii.MCSE_PII[k]),
            file=file,
        )  # fmt: skip
    print(
        "  PII = log10(CI UB/CI LB); NaN when CI LB <= 0 (read the range fraction). "
        "ratio near 1 or atLB/atUB above 0.05: the bounds set the interval.\n",
        file=file,
    )


def _number(x):
    """Full precision; NaN and Inf as MATLAB reads them back."""
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return "%.10g" % x


_PII_COLUMNS = ["estimate", "median", "CI_LB", "CI_UB", "PII", "PII_rangeFrac", "PII_scale", "PII_box", "PII_ratio",
                "fracLB", "fracUB", "MCSE_CI_LB", "MCSE_CI_UB", "MCSE_PII", "bound_LB", "bound_UB"]  # fmt: skip


def subepidemicWritePII(filename, pii, ranks=None):
    """Save the PII metrics as a csv file (``subepidemicWritePII.m``).

    ``subepidemicWritePII(filename, pii)`` writes one Struct, one row per
    parameter. ``subepidemicWritePII(filename, piiList, ranks)`` writes
    several ranked models into one file; the rows then start with the columns
    ``ranked_model`` and ``n_subepidemics``.
    """
    if isinstance(pii, dict):
        pii = [pii]
    withRank = ranks is not None and np.size(ranks) > 0
    if withRank:
        ranks = [int(r) for r in np.asarray(ranks).ravel()]
        if len(ranks) != len(pii):
            raise SubEpiPredictError("SubEpiPredict:PIIRanks", "ranks must have one entry per PII struct.")
    try:
        fid = open(filename, "w", newline="")
    except OSError:
        raise SubEpiPredictError(
            "SubEpiPredict:PIIFile", "The file %s could not be opened for writing." % filename
        ) from None
    with fid:
        if withRank:
            fid.write("ranked_model,n_subepidemics,")
        fid.write("parameter,subepidemic")
        for c in _PII_COLUMNS:
            fid.write(",%s" % c)
        fid.write(",fixed\n")
        for m, S in enumerate(pii):
            subs = np.asarray(S["subepidemic"], dtype=float)
            nSub = int(np.nanmax(subs))
            for k in range(len(S["name"])):
                if withRank:
                    fid.write("%d,%d," % (ranks[m], nSub))
                fid.write("%s,%s" % (S["name"][k], _number(subs[k])))
                for c in _PII_COLUMNS:
                    fid.write(",%s" % _number(S[c][k]))
                fid.write(",%d\n" % int(bool(S["fixed"][k])))
