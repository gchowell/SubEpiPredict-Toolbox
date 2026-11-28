"""The n-sub-epidemic model: kernels, initial state, incidence, event-driven simulation."""

import math

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError

T = np.arange(0, 61, dtype=float)


def sim(flag1, r, p, a, K, thr=0.0, I0=1.0, n=None, onset_fixed=0, t=T):
    r = np.atleast_1d(np.asarray(r, dtype=float))
    n = r.size if n is None else n
    full = lambda v: np.full(n, v, dtype=float) if np.ndim(v) == 0 else np.asarray(v, dtype=float)  # noqa: E731
    return sp.simulateSubepidemic(t, I0, n, onset_fixed, thr, flag1, full(r), full(p), full(a), full(K))


# ---- growth kernels ---------------------------------------------------------


def test_kernels_have_the_six_expressions():
    x, r, p, a, K = 37.0, 0.4, 0.8, 1.7, 500.0
    expected = {
        0: r * x**p,
        1: r * x**p * (1 - x / K),
        2: r * x**p * (1 - x / K) ** a,
        3: r * x * (1 - x / K),
        4: r * x * (1 - (x / K) ** a),
        5: r * x * math.log(K / x),
    }
    for flag1, value in expected.items():
        dx = sp.modifiedLogisticGrowthPatch(0, [x], [r], [p], [a], [K], 1, 0, flag1, [1])
        assert dx[0] == pytest.approx(value, rel=1e-14)


def test_dormant_patch_has_zero_derivative_even_with_nan_parameters():
    dx = sp.modifiedLogisticGrowthPatch(0, [5.0, 1.0], [0.3, np.nan], [0.9, np.nan], [1, 1], [100.0, np.nan], 2, 0, 1, [1, 0])
    assert dx[0] > 0 and dx[1] == 0


def test_kernel_needs_an_activation_mask():
    with pytest.raises(SubEpiPredictError) as e:
        sp.modifiedLogisticGrowthPatch(0, [1.0], [0.1], [1], [1], [10], 1, 0, 1)
    assert identifier_of(e) == "SubEpiPredict:MissingActivationMask"
    with pytest.raises(SubEpiPredictError) as e:
        sp.modifiedLogisticGrowthPatch(0, [1.0], [0.1], [1], [1], [10], 1, 0, 1, [2])
    assert identifier_of(e) == "SubEpiPredict:InvalidActivationMask"
    with pytest.raises(SubEpiPredictError) as e:
        sp.modifiedLogisticGrowthPatch(0, [1.0], [0.1], [1], [1], [10], 1, 0, 7, [1])
    assert identifier_of(e) == "SubEpiPredict:InvalidModelFlag"


# ---- initial state, parameters, incidence ------------------------------------


def test_initial_state_asynchronous_and_synchronous():
    IC, inv, tinv, cinv = sp.initializeSubepidemicState(12, 3, 0)
    assert IC.tolist() == [12, 1, 1] and inv.tolist() == [1, 0, 0]
    assert tinv.tolist() == [0, 0, 0] and cinv.tolist() == [0, 0, 0]
    IC, inv, _, _ = sp.initializeSubepidemicState(12, 3, 1)
    assert IC.tolist() == [4, 4, 4] and inv.tolist() == [1, 1, 1]


@pytest.mark.parametrize(
    "args, identifier",
    [
        ((-1, 2, 0), "SubEpiPredict:InvalidI0"),
        ((np.inf, 2, 0), "SubEpiPredict:InvalidI0"),
        ((1, 0, 0), "SubEpiPredict:InvalidPatchCount"),
        ((1, 2.5, 0), "SubEpiPredict:InvalidPatchCount"),
        ((1, 2, 2), "SubEpiPredict:InvalidOnsetMode"),
    ],
)
def test_initial_state_checks_its_inputs(args, identifier):
    with pytest.raises(SubEpiPredictError) as e:
        sp.initializeSubepidemicState(*args)
    assert identifier_of(e) == identifier


def test_unpack_parameters_layout():
    z = np.arange(1, 11, dtype=float)  # n = 2
    rs, ps, as_, Ks, alpha, d = sp.unpackSubepidemicParameters(z, 2)
    assert rs.tolist() == [1, 2] and ps.tolist() == [3, 4] and as_.tolist() == [5, 6] and Ks.tolist() == [7, 8]
    assert (alpha, d) == (9, 10)
    # a column is accepted
    assert sp.unpackSubepidemicParameters(z.reshape(-1, 1), 2)[3].tolist() == [7, 8]
    with pytest.raises(SubEpiPredictError) as e:
        sp.unpackSubepidemicParameters(z[:-1], 2)
    assert identifier_of(e) == "SubEpiPredict:ParameterDimensionMismatch"


def test_incidence_corrects_only_the_first_row_of_asynchronous_seeds():
    x = np.array([[5.0, 1.0, 1.0], [9.0, 1.0, 1.0], [15.0, 4.0, 1.0]])
    totinc, patch = sp.subepidemicIncidence(x, 0)
    assert patch.tolist() == [[5, 0, 0], [4, 0, 0], [6, 3, 0]]
    assert totinc.tolist() == [5, 4, 9]
    assert np.allclose(patch.sum(axis=1), totinc)
    totinc, patch = sp.subepidemicIncidence(x, 1)
    assert patch[0].tolist() == [5, 1, 1] and totinc[0] == 7


def test_ode_options_are_the_toolbox_tolerances():
    o = sp.subepidemicOdeOptions()
    assert o.RelTol == 1e-8 and o.AbsTol == 1e-10


def test_activation_event_function():
    value, isterminal, direction = sp.subepidemicActivationEvents(0, [30.0, 1.0, 1.0], 25.0, [2, 3])
    assert value.tolist() == [5.0, -24.0] and isterminal.tolist() == [1, 1] and direction.tolist() == [1, 1]


# ---- one sub-epidemic against closed forms -----------------------------------


def test_logistic_closed_form():
    r, K, C0 = 0.25, 5000.0, 3.0
    _, x, totinc, patch, act = sim(3, r, 1, 1, K, I0=C0)
    exact = K / (1 + (K / C0 - 1) * np.exp(-r * T))
    assert np.max(np.abs(x[:, 0] - exact) / exact) < 5e-7
    assert totinc[0] == pytest.approx(C0)
    assert np.allclose(np.cumsum(totinc), x[:, 0])
    assert act.active.tolist() == [True] and act.nSegments == 1


def test_generalized_growth_closed_form():
    r, p, C0 = 0.7, 0.6, 2.0
    _, x, _, _, _ = sim(0, r, p, 1, 20, I0=C0)
    exact = (C0 ** (1 - p) + (1 - p) * r * T) ** (1 / (1 - p))
    assert np.max(np.abs(x[:, 0] - exact) / exact) < 5e-7


def test_gompertz_closed_form():
    r, K, C0 = 0.12, 800.0, 4.0
    _, x, _, _, _ = sim(5, r, 1, 1, K, I0=C0)
    exact = K * np.exp(np.log(C0 / K) * np.exp(-r * T))
    assert np.max(np.abs(x[:, 0] - exact) / exact) < 5e-7


def test_richards_closed_form():
    r, a, K, C0 = 0.3, 0.6, 900.0, 2.0
    _, x, _, _, _ = sim(4, r, 1, a, K, I0=C0)
    Q = (K / C0) ** a - 1
    exact = K / (1 + Q * np.exp(-a * r * T)) ** (1 / a)
    assert np.max(np.abs(x[:, 0] - exact) / exact) < 5e-7


def test_generalized_logistic_reduces_to_logistic_at_p_1():
    a1 = sim(1, 0.25, 1.0, 1, 5000.0, I0=3.0)[2]
    a3 = sim(3, 0.25, 1.0, 1, 5000.0, I0=3.0)[2]
    assert np.allclose(a1, a3, rtol=1e-9, atol=1e-9)


def test_generalized_richards_can_saturate():
    # with a < 1 the curve reaches K in finite time; the simulation must go on
    _, x, totinc, _, _ = sim(2, 0.9, 0.9, 0.5, 300.0, I0=5.0, t=np.arange(0, 201.0))
    assert np.all(np.isfinite(x)) and x[-1, 0] == pytest.approx(300.0, rel=1e-5)
    assert abs(totinc[-1]) < 1e-3


# ---- activation ---------------------------------------------------------------


def test_crossing_time_matches_the_logistic_formula():
    r, K, C0, thr = 0.3, 2000.0, 1.0, 150.0
    _, x, _, _, act = sim(3, [r, r], 1, 1, [K, K], thr=thr, I0=C0)
    tau = math.log((K / C0 - 1) / (K / thr - 1)) / r  # time at which the logistic curve reaches thr
    assert act.time[0] == 0 and act.time[1] == pytest.approx(tau, abs=1e-6)
    assert act.predecessorCount[1] == pytest.approx(thr)
    assert math.isnan(act.predecessorCount[0])
    assert act.initialActive.tolist() == [True, False] and act.active.tolist() == [True, True]
    assert act.nSegments == 2
    # the second sub-epidemic is the same logistic curve started from 1 at tau
    later = T > tau
    exact2 = K / (1 + (K - 1) * np.exp(-r * (T[later] - tau)))
    assert np.max(np.abs(x[later, 1] - exact2) / exact2) < 5e-7
    assert np.all(x[~later, 1] == 1.0)


def test_event_times_are_not_inserted_and_states_are_continuous():
    t, x, totinc, patch, _ = sim(1, [0.3, 0.3], 0.95, 1, [3000.0, 1500.0], thr=200.0)
    assert np.array_equal(t, T) and x.shape == (T.size, 2)
    assert np.all(np.diff(x, axis=0) >= -1e-9)
    assert np.allclose(patch.sum(axis=1), totinc)
    assert totinc[0] == pytest.approx(1.0)  # seeds are not counted as cases


def test_threshold_never_reached_leaves_the_patch_dormant():
    _, x, totinc, patch, act = sim(1, [0.3, 0.3], 0.9, 1, [100.0, 500.0], thr=150.0)
    assert act.active.tolist() == [True, False] and math.isnan(act.time[1])
    assert np.all(x[:, 1] == 1.0) and np.all(patch[:, 1] == 0)
    assert act.nSegments == 1


def test_zero_threshold_starts_all_patches_but_keeps_asynchronous_seeds():
    _, x, totinc, _, act = sim(1, [0.2, 0.3, 0.4], 0.9, 1, [900.0, 800.0, 700.0], thr=0.0, I0=6.0)
    assert x[0].tolist() == [6.0, 1.0, 1.0]
    assert act.active.all() and np.all(act.time == 0) and act.nSegments == 1
    assert totinc[0] == pytest.approx(6.0)


def test_synchronous_start_splits_the_initial_value():
    _, x, totinc, _, act = sim(1, [0.2, 0.3, 0.4], 0.9, 1, [900.0, 800.0, 700.0], thr=123.0, I0=6.0, onset_fixed=1)
    assert x[0].tolist() == [2.0, 2.0, 2.0] and totinc[0] == pytest.approx(6.0)
    assert act.initialActive.all() and act.nSegments == 1


def test_literal_threshold_rule_for_three_patches():
    # threshold <= 1: the third patch starts at time zero because the dormant
    # second patch holds its seed of 1; the second waits for the first.
    _, x, _, _, act = sim(1, [0.2, 0.2, 0.2], 1.0, 1, [500.0, 500.0, 500.0], thr=0.8, I0=0.5)
    assert act.time[2] == 0 and act.time[1] > 0
    # threshold above 1: a cascade
    _, _, _, _, act = sim(1, [0.3, 0.3, 0.3], 1.0, 1, [500.0, 500.0, 500.0], thr=50.0, I0=1.0)
    assert 0 < act.time[1] < act.time[2] and act.nSegments == 3


def test_two_point_grid_and_single_point():
    t, x, totinc, _, _ = sim(3, 0.5, 1, 1, 100.0, I0=2.0, t=np.array([0.0, 1.0]))
    assert x.shape == (2, 1) and totinc.size == 2
    t, x, totinc, _, act = sim(3, 0.5, 1, 1, 100.0, I0=2.0, t=np.array([0.0]))
    assert x.tolist() == [[2.0]] and act.nSegments == 0


def test_results_do_not_depend_on_the_grid():
    fine = sim(1, [0.3, 0.3], 0.95, 1, [3000.0, 1500.0], thr=200.0, t=np.arange(0, 60.5, 0.5))[1]
    coarse = sim(1, [0.3, 0.3], 0.95, 1, [3000.0, 1500.0], thr=200.0, t=np.arange(0, 61.0))[1]
    assert np.max(np.abs(fine[::2] - coarse) / coarse) < 1e-6


@pytest.mark.parametrize(
    "change, identifier",
    [
        (dict(t=np.array([0.0, 1.0, 1.0])), "SubEpiPredict:InvalidTimeGrid"),
        (dict(t=np.array([0.0, np.nan])), "SubEpiPredict:InvalidTimeGrid"),
        (dict(thr=-1.0), "SubEpiPredict:InvalidThreshold"),
        (dict(flag1=6), "SubEpiPredict:InvalidModelFlag"),
        (dict(rs=[0.1]), "SubEpiPredict:ParameterDimensionMismatch"),
        (dict(rs=[0.1, np.inf]), "SubEpiPredict:InvalidGrowthParameters"),
        (dict(ode={"Events": 1}), "SubEpiPredict:ReservedEventsOption"),
        (dict(ode={"Jacobian": 1}), "SubEpiPredict:UnsupportedOdeOption"),
        (dict(ode=3), "SubEpiPredict:InvalidOdeOptions"),
    ],
)
def test_simulation_checks_its_inputs(change, identifier):
    a = dict(t=T, thr=10.0, flag1=1, rs=[0.2, 0.2], ode=None)
    a.update(change)
    with pytest.raises(SubEpiPredictError) as e:
        sp.simulateSubepidemic(a["t"], 1, 2, 0, a["thr"], a["flag1"], a["rs"], [1, 1], [1, 1], [100, 100], a["ode"])
    assert identifier_of(e) == identifier


def test_looser_tolerance_changes_the_last_digits_only():
    tight = sim(1, 0.5, 0.8, 1, 8000.0)[1]
    loose = sp.simulateSubepidemic(T, 1, 1, 0, 0, 1, [0.5], [0.8], [1], [8000.0], sp.Struct(RelTol=1e-4, AbsTol=1e-6))[1]
    assert 0 < np.max(np.abs(tight - loose) / tight) < 1e-2


def test_simulate_refit_reports_a_failed_curve_instead_of_raising():
    ok, x, totinc, patch, message = sp.subepidemicSimulateRefit(T, 1, 1, 0, 0, 1, [0.5], [0.8], [1], [8000.0])
    assert ok and message == "" and totinc.size == T.size
    # super-exponential growth blows up in finite time
    ok, x, totinc, patch, message = sp.subepidemicSimulateRefit(np.arange(0, 400.0), 5, 1, 0, 0, 0, [5.0], [3.0], [1], [20.0])
    assert not ok and message
    # wrong sizes are still an error
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicSimulateRefit(T, 1, 2, 0, 0, 1, [0.5], [0.8], [1], [8000.0])


def coupled_event_solution(t, I0, thr, flag1, rs, ps, as_, Ks):
    """The model as MATLAB integrates it: all sub-epidemics together, stopping
    at each upward threshold crossing and restarting with the new mask."""
    from scipy.integrate import solve_ivp

    n = len(rs)
    state = np.ones(n)
    state[0] = I0
    active = np.zeros(n, dtype=bool)
    active[0] = True
    for j in range(1, n):
        active[j] = state[j - 1] >= thr
    x = np.full((t.size, n), np.nan)
    x[0] = state
    now = t[0]
    while now < t[-1]:
        mask = active.copy()
        watched = [j for j in range(1, n) if not mask[j]]
        events = []
        for w in watched:
            def event(tt, xx, w=w):
                return xx[w - 1] - thr

            event.terminal = True
            event.direction = 1
            events.append(event)
        sol = solve_ivp(
            lambda tt, xx: sp.modifiedLogisticGrowthPatch(tt, xx, rs, ps, as_, Ks, n, thr, flag1, mask),
            (now, t[-1]), state, method="LSODA", rtol=1e-11, atol=1e-12, dense_output=True, events=events or None,
        )  # fmt: skip
        rows = (t > now) & (t <= sol.t[-1])
        x[rows] = sol.sol(t[rows]).T
        now = sol.t[-1]
        state = sol.y[:, -1].copy()
        for k, w in enumerate(watched):
            if len(sol.t_events[k]):
                active[w] = True
        for j in range(1, n):
            if not active[j] and state[j - 1] >= thr:
                active[j] = True
    return x


@pytest.mark.parametrize("flag1, thr", [(1, 100.0), (1, 3000.0), (1, 5000.0), (2, 400.0), (4, 300.0), (5, 250.0), (1, 0.5)])
def test_agrees_with_an_event_driven_solve_of_the_coupled_system(flag1, thr):
    t = np.arange(0, 221, dtype=float)
    rs, ps, as_, Ks = [0.18, 0.18, 0.18], [0.98, 0.98, 0.98], [1.0, 1.3, 0.8], [10000.0, 5000.0, 1000.0]
    x = sp.simulateSubepidemic(t, 1, 3, 0, thr, flag1, rs, ps, as_, Ks)[1]
    reference = coupled_event_solution(t, 1.0, thr, flag1, rs, ps, as_, Ks)
    assert np.max(np.abs(x - reference) / np.maximum(1, np.abs(reference))) < 1e-6
