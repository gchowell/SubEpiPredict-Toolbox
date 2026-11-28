"""Search space, unit map, start points and the order of interchangeable labels."""

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError

DATA = np.array([4.0, 6.0, 9.0, 15.0, 30.0, 50.0, 41.0, 22.0])  # sum 177, max 50


def test_legacy_bounds_of_the_generalized_logistic_model():
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(2, DATA, 1, 3)
    # layout [r1 r2 p1 p2 a1 a2 K1 K2 alpha d]
    assert LB.tolist() == [0.025, 0.025, 0, 0, 1, 1, 20, 20, 1e-8, 0]
    assert UB.tolist() == [500, 500, 1, 1, 1, 1, 1e10, 1e10, 1e5, 0]  # r: max(50*6, 10*50)
    assert isLog.tolist() == [True, True, False, False, False, False, True, True, True, False]
    assert z0.tolist() == [0.1, 0.1, 0.9, 0.9, 1, 1, 177, 177, 1, 0]


@pytest.mark.parametrize(
    "flag1, p, a, K",
    [(0, (0, 1), (1, 1), (20, 20)), (1, (0, 1), (1, 1), (20, 1e10)), (2, (0, 1), (0, 10), (20, 1e10)),
     (3, (1, 1), (1, 1), (20, 1e10)), (4, (1, 1), (0, 10), (20, 1e10)), (5, (1, 1), (1, 1), (20, 1e10))],
)
def test_bounds_of_every_growth_model(flag1, p, a, K):
    LB, UB, _, z0 = sp.subepidemicSearchSpace(1, DATA, flag1, 0)
    assert (LB[1], UB[1]) == p and (LB[2], UB[2]) == a and (LB[3], UB[3]) == K
    assert np.all(z0 >= LB) and np.all(z0 <= UB)
    assert LB[4:].tolist() == [0, 0] and UB[4:].tolist() == [0, 0]  # least squares: alpha and d fixed


@pytest.mark.parametrize("method1, lb, ub", [(1, [0, 0], [0, 0]), (4, [1e-8, 0], [1e5, 0]), (5, [1e-8, 0.2], [1e5, 100]),
                                             (6, [0, 0], [0, 0])])  # fmt: skip
def test_bounds_of_the_error_parameters(method1, lb, ub):
    LB, UB, _, _ = sp.subepidemicSearchSpace(1, DATA, 1, method1)
    assert LB[4:].tolist() == lb and UB[4:].tolist() == ub


def test_settings_replace_bounds_and_scales():
    S = sp.subepidemicFitSettings(dict(rUpperBoundFactor=0, KUB=[2e6, 5e6, 9e9], pLB=1, pUB=1, logScaleK=False))
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(2, DATA, 1, 0, S)
    assert UB[0] == 300  # the legacy bound alone: 50*max(first two observations)
    assert UB[6:8].tolist() == [2e6, 5e6]
    assert LB[2:4].tolist() == [1, 1] and UB[2:4].tolist() == [1, 1] and z0[2:4].tolist() == [1, 1]
    assert not isLog[6] and isLog[0]
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicSearchSpace(3, DATA, 1, 0, sp.subepidemicFitSettings(dict(KUB=[2e6, 5e6])))
    assert identifier_of(e) == "SubEpiPredict:FitSettings"
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicSearchSpace(1, DATA, 1, 0, sp.subepidemicFitSettings(dict(rUB=1e-6)))
    assert identifier_of(e) == "SubEpiPredict:ReversedBounds"


def test_search_space_checks_its_inputs():
    for args, identifier in (
        ((0, DATA, 1, 0), "SubEpiPredict:InvalidPatchCount"),
        ((1, [1.0], 1, 0), "SubEpiPredict:InvalidData"),
        ((1, DATA, 9, 0), "SubEpiPredict:InvalidModelFlag"),
        ((1, DATA, 1, 7), "SubEpiPredict:InvalidMethod"),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicSearchSpace(*args)
        assert identifier_of(e) == identifier


def make_map(n=2, method1=3):
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(n, DATA, 1, method1)
    return sp.subepidemicUnitMap("make", dict(z0=z0, LB=LB, UB=UB, isLog=isLog))


def test_unit_map_round_trip_and_fixed_entries():
    m = make_map()
    assert m.nz == 10 and m.nFree == 7 and m.free.tolist() == [0, 1, 2, 3, 6, 7, 8]
    z = np.array([3.0, 40.0, 0.2, 0.7, 1, 1, 5e3, 8e7, 12.0, 0])
    u = sp.subepidemicUnitMap("toU", z, m)
    assert u.shape == (1, 7) and np.all((u >= 0) & (u <= 1))
    assert np.allclose(sp.subepidemicUnitMap("toZ", u[0], m), z, rtol=1e-12)
    # log coordinate: the middle of the unit interval is the geometric mean of the bounds
    mid = sp.subepidemicUnitMap("toZ", np.full(7, 0.5), m)
    assert mid[6] == pytest.approx(np.sqrt(20 * 1e10)) and mid[2] == pytest.approx(0.5)
    # corners are the bounds; values outside the box are clipped
    low = sp.subepidemicUnitMap("toZ", np.zeros(7), m)
    high = sp.subepidemicUnitMap("toZ", np.full(7, 2.0), m)
    assert np.allclose(low, m.LB, rtol=1e-13) and np.all(low >= m.LB)
    assert np.allclose(high, m.UB, rtol=1e-13) and np.all(high <= m.UB)
    assert np.all(sp.subepidemicUnitMap("toU", m.UB * 10, m) == 1)


def test_unit_map_errors():
    m = make_map()
    for action, arg, identifier in (
        ("toZ", np.zeros(3), "SubEpiPredict:ParameterDimensionMismatch"),
        ("toU", np.zeros((2, 4)), "SubEpiPredict:ParameterDimensionMismatch"),
        ("other", None, "SubEpiPredict:UnitMapAction"),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicUnitMap(action, arg, m)
        assert identifier_of(e) == identifier
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicUnitMap("make", dict(z0=[1, 1], LB=[2, 0], UB=[1, 1], isLog=[0, 0]))
    assert identifier_of(e) == "SubEpiPredict:ReversedBounds"


def test_start_points_cover_the_box_and_end_with_the_seed():
    m = make_map()
    sp.rng(3)
    starts = sp.subepidemicStartPoints(m, 20, 3)
    assert starts.shape == (24, 10)
    assert np.array_equal(starts[-1], m.z0)
    assert np.all(starts >= m.LB) and np.all(starts <= m.UB)
    fixed = m.UB <= m.LB
    assert np.all(starts[:, fixed] == m.z0[fixed])
    # Latin hypercube in unit coordinates: one point per stratum in every free coordinate
    U = sp.subepidemicUnitMap("toU", starts[:20], m)
    for j in range(m.nFree):
        assert sorted(np.floor(U[:, j] * 20).astype(int).tolist()) == list(range(20))
    # K is sampled log-uniformly: about half of the points are below the geometric mean
    assert 5 <= np.sum(starts[:20, 6] < np.sqrt(20 * 1e10)) <= 15
    # same seed, same points
    sp.rng(3)
    assert np.array_equal(sp.subepidemicStartPoints(m, 20, 3), starts)


def test_start_points_special_cases():
    m = make_map()
    assert sp.subepidemicStartPoints(m, 0, 0).tolist() == [m.z0.tolist()]
    seed = m.z0.copy()
    seed[0] = 7.0
    starts = sp.subepidemicStartPoints(m, 2, 2, seed)
    assert starts[-1, 0] == 7.0 and starts.shape == (5, 10)
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicStartPoints(m, 1, 0, [1, 2, 3])
    assert identifier_of(e) == "SubEpiPredict:InvalidSeed"


def test_lhsdesign_is_a_latin_hypercube():
    X = sp.lhsdesign(7, 3, np.random.default_rng(0))
    assert X.shape == (7, 3)
    for j in range(3):
        assert sorted(np.floor(X[:, j] * 7).astype(int).tolist()) == list(range(7))


# ---- order of interchangeable labels -----------------------------------------


def P(r, K, n):
    return np.concatenate([r, np.arange(n) + 0.1, np.arange(n) + 0.5, K, [2.0, 0.0]])


def test_synchronous_models_are_sorted_by_r_then_K():
    p = P([0.5, 0.2, 0.2], [100.0, 300.0, 200.0], 3)
    out, group, changed = sp.subepidemicOrderLabels(p, 3, 1, 0, 5)
    assert group.tolist() == [0, 1, 2] and changed.tolist() == [True]
    assert out[0:3].tolist() == [0.2, 0.2, 0.5]
    assert out[9:12].tolist() == [200.0, 300.0, 100.0]
    assert out[3:6].tolist() == [2.1, 1.1, 0.1]  # p moves with its sub-epidemic
    assert out[-2:].tolist() == [2.0, 0.0]


def test_asynchronous_groups():
    two = P([0.5, 0.2], [100.0, 300.0], 2)
    # threshold met at time zero and I0 = 1: both interchangeable
    out, group, changed = sp.subepidemicOrderLabels(two, 2, 0, 0, 1)
    assert group.tolist() == [0, 1] and out[0:2].tolist() == [0.2, 0.5]
    # I0 > 1 with two sub-epidemics: they start from different values, nothing to do
    out, group, changed = sp.subepidemicOrderLabels(two, 2, 0, 0, 8)
    assert group.size == 0 and not changed[0] and np.array_equal(out, two)
    # threshold not met at time zero: the onsets set the order
    out, group, changed = sp.subepidemicOrderLabels(two, 2, 0, 50, 1)
    assert group.size == 0 and np.array_equal(out, two)
    # three sub-epidemics, I0 > 1, threshold <= 1: sub-epidemics 2 and 3 are interchangeable
    three = P([0.9, 0.5, 0.2], [100.0, 300.0, 200.0], 3)
    out, group, changed = sp.subepidemicOrderLabels(three, 3, 0, 0.5, 8)
    assert group.tolist() == [1, 2] and out[0:3].tolist() == [0.9, 0.2, 0.5]


def test_order_labels_keeps_the_curve_and_handles_matrices():
    t = np.arange(0, 40.0)
    p = P([0.4, 0.25], [900.0, 400.0], 2)
    p[2:6] = [0.9, 0.8, 1, 1]
    q, _, _ = sp.subepidemicOrderLabels(p, 2, 1, 0, 6)
    a = sp.simulateSubepidemic(t, 6, 2, 1, 0, 1, p[0:2], p[2:4], p[4:6], p[6:8])[2]
    b = sp.simulateSubepidemic(t, 6, 2, 1, 0, 1, q[0:2], q[2:4], q[4:6], q[6:8])[2]
    assert np.allclose(a, b, rtol=1e-12)
    rows = np.vstack([p, q, np.full(10, np.nan)])
    out, _, changed = sp.subepidemicOrderLabels(rows, 2, 1, 0, 6)
    assert changed.tolist() == [True, False, False] and np.array_equal(out[0], out[1])
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicOrderLabels(np.zeros(7), 2, 1, 0, 6)
