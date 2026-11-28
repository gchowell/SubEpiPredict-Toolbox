"""Search space of one candidate model: bounds, coordinate scales, start points.

Ports of ``subepidemicSearchSpace.m``, ``subepidemicUnitMap.m``,
``subepidemicStartPoints.m`` and ``subepidemicOrderLabels.m``.

Layout of every parameter vector: ``[r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]``,
``4*n + 2`` entries.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError
from .settings import resolve_fit_settings
from .utils import Struct, get_rng

_REALMIN = np.finfo(float).tiny


def _check_patch_count(npatches):
    ok = isinstance(npatches, (int, float, np.integer, np.floating)) and not isinstance(npatches, (bool, np.bool_))
    if not (ok and math.isfinite(npatches) and npatches >= 1 and npatches == int(npatches)):
        raise SubEpiPredictError("SubEpiPredict:InvalidPatchCount", "npatches must be a positive integer.")
    return int(npatches)


def _localBound(settings, name, legacy, npatches):
    """One bound per sub-epidemic: the legacy value, or the setting when given."""
    b = np.zeros(npatches) + legacy
    v = settings.get(name)
    if v is not None and np.size(v) > 0:
        v = np.asarray(v, dtype=float).ravel()
        if v.size == 1:
            b[:] = v[0]
        elif v.size >= npatches:
            b = v[:npatches].copy()
        else:
            raise SubEpiPredictError(
                "SubEpiPredict:FitSettings",
                "Setting %s has %d entries; %d sub-epidemics need a scalar or at least %d values."
                % (name, v.size, npatches, npatches),
            )
    return b


def subepidemicSearchSpace(npatches, data, flag1, method1, settings=None):
    """Bounds, coordinate scales and seed of one candidate.

    Port of ``subepidemicSearchSpace.m``.

    Returns ``(LB, UB, isLog, z0)``.

    ``LB``, ``UB``  bounds. r: from ``mean(|first two observations|)/200`` to
        the larger of ``50*max(|first two observations|)`` and
        ``settings.rUpperBoundFactor`` (10) times the largest count, for the
        kernels with p free (GGM, GLM, GRM). With ``settings.rBoundsByKernel``
        (default True) the kernels with p = 1 use ``settings.rRateBounds``
        (default 1e-3 to 10 per time step) instead: logistic and Gompertz
        both ends, Richards the lower end only. p: 0 to 1
        (GGM, GLM, GRM), fixed at 1 otherwise. a: 0 to 10 for GRM and
        Richards, fixed at 1 otherwise (also for Gompertz, whose equation does
        not contain a). K: 20 to 1e10 (fixed at 20 for GGM, which has no K).
        alpha: 1e-8 to 1e5 for methods 3, 4, 5, fixed at 0 otherwise. d: 0.2
        to 100 for method 5, fixed at 0 otherwise. A parameter with
        ``LB == UB`` is fixed. The settings ``rLB`` ... ``KUB`` replace a
        bound. Reversed bounds are an error.
    ``isLog``  True for the coordinates searched on a log scale: r, K and
        alpha, when the setting is on, the lower bound is positive and the
        parameter is free.
    ``z0``  the legacy seed (r = 0.1, p = 0.9 or 1, a = 1, K = sum of the
        data, alpha = 1, d = 1), projected into the bounds.

    ``data`` is the calibration incidence (the unsmoothed series).
    """
    settings = resolve_fit_settings(settings)
    npatches = _check_patch_count(npatches)
    try:
        data = np.asarray(data, dtype=float)
        ok = data.size >= 2 and (data.ndim <= 1 or (data.ndim == 2 and 1 in data.shape))
        ok = ok and bool(np.all(np.isfinite(data)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidData", "data must be a finite real vector with at least two observations."
        )
    data = data.ravel()
    if not (np.ndim(flag1) == 0 and flag1 in (0, 1, 2, 3, 4, 5)):
        raise SubEpiPredictError("SubEpiPredict:InvalidModelFlag", "flag1 must be an integer from 0 to 5.")
    if not (np.ndim(method1) == 0 and method1 in (0, 1, 2, 3, 4, 5, 6)):
        raise SubEpiPredictError("SubEpiPredict:InvalidMethod", "method1 must be an integer from 0 to 6.")

    # ---- legacy bounds of getbounds.m ---------------------------------------
    rlb = float(np.mean(np.abs(data[:2]))) / 200
    rub = float(np.max(np.abs(data[:2]))) * 50
    # The legacy upper bound of r depends on the first two observations only. A
    # multiple of the largest count scales with the data and does not constrain
    # a sensible fit.
    if settings.get("rUpperBoundFactor") is not None:
        rub = max(rub, settings.rUpperBoundFactor * float(np.max(np.abs(data))))
    # Kernels with p fixed at 1 (logistic, Richards, Gompertz): r is a rate per
    # time step, while the bounds above are counts (r has units
    # cases^(1-p)/time only when p is free). For a window that starts with
    # large counts the count-based lower bound excludes every plausible rate
    # (r >= 10 per step for first counts near 2000), and the fit escapes to a
    # curve that saturates below the counts already observed. With
    # rBoundsByKernel these kernels get rRateBounds instead: both ends for the
    # logistic and Gompertz kernels, the lower end only for Richards, whose r
    # trades off against a as a -> 0 and can exceed any per-step rate.
    if settings.rBoundsByKernel and flag1 in (3, 4, 5):
        rlb = float(settings.rRateBounds[0])
        if flag1 != 4:
            rub = float(settings.rRateBounds[1])
    Kmax = 10000000000.0
    if flag1 == 0:  # GGM
        pb, ab, Kb = (0, 1), (1, 1), (20, 20)
    elif flag1 == 1:  # GLM
        pb, ab, Kb = (0, 1), (1, 1), (20, Kmax)
    elif flag1 == 2:  # GRM
        pb, ab, Kb = (0, 1), (0, 10), (20, Kmax)
    elif flag1 == 3:  # Logistic
        pb, ab, Kb = (1, 1), (1, 1), (20, Kmax)
    elif flag1 == 4:  # Richards
        pb, ab, Kb = (1, 1), (0, 10), (20, Kmax)
    else:  # Gompertz: dC/dt = r*C*log(K/C) has no a
        pb, ab, Kb = (1, 1), (1, 1), (20, Kmax)
    rLB = _localBound(settings, "rLB", rlb, npatches)
    rUB = _localBound(settings, "rUB", rub, npatches)
    pLB = _localBound(settings, "pLB", pb[0], npatches)
    pUB = _localBound(settings, "pUB", pb[1], npatches)
    aLB = _localBound(settings, "aLB", ab[0], npatches)
    aUB = _localBound(settings, "aUB", ab[1], npatches)
    KLB = _localBound(settings, "KLB", Kb[0], npatches)
    KUB = _localBound(settings, "KUB", Kb[1], npatches)

    # ---- legacy bounds of the error-structure parameters (alpha, d) ---------
    if method1 in (0, 1, 2, 6):
        LBe, UBe = [0, 0], [0, 0]
    elif method1 in (3, 4):
        LBe, UBe = [1e-8, 0], [1e5, 0]
    else:
        LBe, UBe = [1e-8, 0.2], [1e5, 1e2]

    LB = np.concatenate([rLB, pLB, aLB, KLB, LBe]).astype(float)
    UB = np.concatenate([rUB, pUB, aUB, KUB, UBe]).astype(float)
    if np.any(LB > UB):
        raise SubEpiPredictError(
            "SubEpiPredict:ReversedBounds",
            "A lower bound exceeds its upper bound. Check the bound settings (rLB ... KUB).",
        )

    # ---- coordinate scales --------------------------------------------------
    isLog = np.zeros(4 * npatches + 2, dtype=bool)
    isLog[0:npatches] = settings.logScaleR
    isLog[3 * npatches : 4 * npatches] = settings.logScaleK
    isLog[4 * npatches] = settings.logScaleAlpha
    isLog = isLog & (LB > 0) & (UB > LB)

    # ---- legacy seed --------------------------------------------------------
    I0 = data[0]
    total = float(np.sum(data))
    rs = np.zeros(npatches) + 0.1
    ps = np.zeros(npatches) + 0.9
    as_ = np.ones(npatches)
    Ks = np.zeros(npatches) + total
    if flag1 in (3, 4, 5):  # logistic, Richards and Gompertz have p = 1
        ps[:] = 1
    if flag1 == 5:
        with np.errstate(all="ignore"):
            rs[:] = 1 - np.divide(I0, total)
            as_[:] = np.divide(rs[0], np.log(np.divide(total, I0)))
    z0 = np.concatenate([rs, ps, as_, Ks, [1.0, 1.0]])
    z0 = np.minimum(np.maximum(z0, LB), UB)
    bad = ~np.isfinite(z0)
    z0[bad] = LB[bad]
    return LB, UB, isLog, z0


# <============================================================================>
# <================================= Unit map =================================>
# <============================================================================>


def _map_toZ(u, m):
    """One point in unit coordinates -> the full vector in original units."""
    u = np.minimum(np.maximum(np.asarray(u, dtype=float).ravel(), 0.0), 1.0)
    g = m.gLB[m.free] + u * m.span[m.free]
    lg = m.isLog[m.free]
    if np.any(lg):
        g[lg] = np.exp(g[lg])
    out = m.z0.copy()
    out[m.free] = np.minimum(np.maximum(g, m.LB[m.free]), m.UB[m.free])
    return out


def _map_toU(Z, m):
    """Parameter vectors (rows) -> unit coordinates of the free entries (rows)."""
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(1, -1)
    if Z.shape[1] != m.nz:
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Each row must hold %d parameters; received %d columns." % (m.nz, Z.shape[1]),
        )
    G = Z.copy()
    if np.any(m.isLog):
        G[:, m.isLog] = np.log(np.maximum(Z[:, m.isLog], _REALMIN))
    out = (G[:, m.free] - m.gLB[m.free]) / m.span[m.free]
    return np.minimum(np.maximum(out, 0.0), 1.0)


def subepidemicUnitMap(action, arg, map=None):
    """Map between parameter vectors and unit coordinates.

    Port of ``subepidemicUnitMap.m``. The optimizer sees only the FREE entries
    (``UB > LB``) of ``z = [r, p, a, K, alpha, d]``, each mapped to [0, 1]:
    ``u = (z - LB)/(UB - LB)`` for a linear coordinate and
    ``u = (log z - log LB)/(log UB - log LB)`` for a log coordinate.

    ``subepidemicUnitMap('make', dict(z0=z0, LB=LB, UB=UB, isLog=isLog))``
        builds the map (a Struct). ``z0`` is projected into the bounds; its
        fixed entries are the values every mapped vector keeps. ``map.free``
        holds the ZERO-BASED indices of the free entries.
    ``subepidemicUnitMap('toU', Z, map)``
        ``Z`` holds one parameter vector per row; the result holds the unit
        coordinates of the free entries, one point per row, clipped to [0, 1].
    ``subepidemicUnitMap('toZ', u, map)``
        ``u`` is ONE point in unit coordinates; the result is the full vector
        in original units, inside the bounds.
    """
    if action == "make":
        z0 = np.asarray(arg["z0"], dtype=float).ravel().copy()
        LB = np.asarray(arg["LB"], dtype=float).ravel().copy()
        UB = np.asarray(arg["UB"], dtype=float).ravel().copy()
        isLog = np.asarray(arg["isLog"]).astype(bool).ravel().copy()
        nz = LB.size
        if any(v.size != nz for v in (z0, UB, isLog)):
            raise SubEpiPredictError(
                "SubEpiPredict:ParameterDimensionMismatch", "Seed, bounds and scale flags must have the same length."
            )
        if np.any(~np.isfinite(LB)) or np.any(~np.isfinite(UB)) or np.any(LB > UB):
            raise SubEpiPredictError(
                "SubEpiPredict:ReversedBounds",
                "Bounds must be finite with LB <= UB. Reversed bounds are not swapped.",
            )
        z0 = np.minimum(np.maximum(z0, LB), UB)
        if np.any(~np.isfinite(z0)):
            raise SubEpiPredictError("SubEpiPredict:InvalidSeed", "The seed is not a finite real vector.")
        isLog = isLog & (LB > 0) & (UB > LB)
        gLB = LB.copy()
        gUB = UB.copy()
        gLB[isLog] = np.log(LB[isLog])
        gUB[isLog] = np.log(UB[isLog])
        free = np.flatnonzero(UB > LB)
        return Struct(z0=z0, LB=LB, UB=UB, isLog=isLog, free=free, gLB=gLB, span=gUB - gLB, nz=nz, nFree=free.size)

    if action == "toU":
        return _map_toU(arg, map)

    if action == "toZ":
        u = np.asarray(arg, dtype=float).ravel()
        if u.size != map.nFree:
            raise SubEpiPredictError(
                "SubEpiPredict:ParameterDimensionMismatch",
                "Expected %d unit coordinates; received %d." % (map.nFree, u.size),
            )
        return _map_toZ(u, map)

    raise SubEpiPredictError("SubEpiPredict:UnitMapAction", "Action must be 'make', 'toU' or 'toZ'.")


# <============================================================================>
# <================================ Start points ==============================>
# <============================================================================>


def lhsdesign(n, p, rng=None, iterations=50):
    """Latin hypercube sample, best of ``iterations`` by the maximin criterion.

    Stands in for ``lhsdesign(n,p,'criterion','maximin','iterations',50)``:
    each column is a random permutation of the n strata with a uniform draw
    inside each stratum; among the candidate designs the one with the largest
    minimum distance between points is kept.
    """
    rng = get_rng(rng)
    best = None
    best_score = -math.inf
    for _ in range(max(1, int(iterations))):
        X = np.empty((n, p))
        for j in range(p):
            X[:, j] = (rng.permutation(n) + rng.random(n)) / n
        if n < 2:
            return X
        diff = X[:, None, :] - X[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", diff, diff)
        score = d2[np.triu_indices(n, 1)].min()
        if score > best_score:
            best_score = score
            best = X
    return best


def unique_rows_stable(A):
    """``unique(A, 'rows', 'stable')``: distinct rows in order of first appearance.

    Returns ``(rows, index)``.
    """
    A = np.asarray(A, dtype=float)
    if A.ndim == 1:
        A = A.reshape(1, -1)
    keep = []
    seen = set()
    for i in range(A.shape[0]):
        key = (A[i] + 0.0).tobytes()  # + 0.0 turns -0.0 into 0.0
        if key not in seen:
            seen.add(key)
            keep.append(i)
    keep = np.array(keep, dtype=int)
    return A[keep].reshape(-1, A.shape[1]), keep


def subepidemicStartPoints(map, nLHS, nJitter, seed=None, rng=None):
    """Start points for the multi-start search, one per row, original units.

    Port of ``subepidemicStartPoints.m``.

    ``map``      from ``subepidemicUnitMap('make', ...)``
    ``nLHS``     Latin-hypercube points over the whole box, drawn in UNIT
                 coordinates, so a log-scale parameter is sampled log-uniformly
    ``nJitter``  copies of ``seed`` with every free parameter multiplied by
                 ``10**N(0, 0.1)``; a free parameter that is zero gets an
                 additive perturbation of 1% of its range instead
    ``seed``     a parameter vector (default ``map.z0``); it is always the
                 last row. Fixed parameters keep the value in ``map.z0``.
    ``rng``      random stream (default: the toolbox stream)

    Duplicate rows are removed.
    """
    rng = get_rng(rng)
    nLHS = int(nLHS)
    nJitter = int(nJitter)
    if seed is None or np.size(seed) == 0:
        seed = map.z0
    try:
        seed = np.asarray(seed, dtype=float).ravel().copy()
        ok = seed.size == map.nz and bool(np.all(np.isfinite(seed)))
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidSeed", "The seed must be a finite real vector with %d entries." % map.nz
        )
    seed = np.minimum(np.maximum(seed, map.LB), map.UB)
    fixed = map.UB <= map.LB
    seed[fixed] = map.z0[fixed]

    # Latin hypercube over the unit box of the free coordinates
    S1 = np.zeros((nLHS, map.nz))
    if nLHS > 0 and map.nFree > 0:
        X = lhsdesign(nLHS, map.nFree, rng)
        for i in range(nLHS):
            S1[i, :] = _map_toZ(X[i, :], map)
    elif nLHS > 0:
        S1 = np.tile(map.z0, (nLHS, 1))

    # jittered copies of the seed
    S2 = np.zeros((nJitter, map.nz))
    if nJitter > 0:
        span = map.UB - map.LB
        jitter = 10.0 ** (0.1 * rng.standard_normal((nJitter, map.nz)))
        S2 = np.tile(seed, (nJitter, 1)) * jitter
        zeroCols = seed == 0
        if np.any(zeroCols):
            S2[:, zeroCols] = seed[zeroCols] + 0.01 * rng.standard_normal((nJitter, int(zeroCols.sum()))) * span[zeroCols]
        S2 = np.minimum(np.maximum(S2, map.LB), map.UB)
        S2[:, fixed] = map.z0[fixed]

    starts = np.vstack([S1, S2, seed.reshape(1, -1)])
    starts, _ = unique_rows_stable(starts)
    keep = np.all(np.isfinite(starts), axis=1)
    starts = starts[keep, :]
    if starts.shape[0] == 0:
        starts = seed.reshape(1, -1)
    return starts


# <============================================================================>
# <============================== Order of labels =============================>
# <============================================================================>


def subepidemicOrderLabels(P, npatches, onset_fixed, onset_thr, I0):
    """Fixed order for sub-epidemics that start together.

    Port of ``subepidemicOrderLabels.m``. Sub-epidemics that start at time
    zero from the same value are interchangeable: swapping their parameters
    gives the same curve and likelihood. Within such a group they are put in
    increasing order of r (ties: increasing K).

    Labels are interchangeable for all sub-epidemics when ``onset_fixed = 1``;
    when ``onset_fixed = 0`` and the onset threshold is already met at time
    zero (``onset_thr <= I0`` for 2 sub-epidemics, ``<= min(I0, 1)`` for
    more): all sub-epidemics if ``I0 = 1``, otherwise sub-epidemics 2 to n
    (3 or more sub-epidemics only).

    ``P`` holds one parameter vector per row (a single vector is accepted and
    returned as a vector). Returns ``(P, group, changed)``: ``group`` the
    ZERO-BASED indices of the interchangeable sub-epidemics (empty: none),
    ``changed`` one flag per row.
    """
    npatches = int(npatches)
    nz = 4 * npatches + 2
    P = np.array(P, dtype=float, copy=True)
    wasVector = P.ndim == 1
    if wasVector:
        P = P.reshape(1, -1)
    if P.ndim != 2 or P.shape[1] != nz:
        raise SubEpiPredictError(
            "SubEpiPredict:ParameterDimensionMismatch",
            "Each row must hold %d parameters; received %d columns." % (nz, P.shape[-1]),
        )
    changed = np.zeros(P.shape[0], dtype=bool)
    group = np.zeros(0, dtype=int)
    if npatches >= 2:
        if onset_fixed == 1:
            group = np.arange(npatches)
        else:
            levelAtStart = I0 if npatches == 2 else min(I0, 1)
            if onset_thr <= levelAtStart:
                if I0 == 1:
                    group = np.arange(npatches)
                elif npatches >= 3:
                    group = np.arange(1, npatches)
    if group.size >= 2:
        rCols = group
        KCols = 3 * npatches + group
        for i in range(P.shape[0]):
            if np.any(~np.isfinite(P[i, np.concatenate([rCols, KCols])])):
                continue
            order = np.lexsort((P[i, KCols], P[i, rCols]))  # by r, ties by K
            if not np.array_equal(order, np.arange(group.size)):
                for block in range(4):
                    P[i, block * npatches + group] = P[i, block * npatches + group[order]]
                changed[i] = True
    if wasVector:
        P = P.ravel()
    return P, group, changed
