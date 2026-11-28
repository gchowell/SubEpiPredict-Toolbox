"""The multi-start fit of one candidate, the screen and the bootstrap."""

import math
import os

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError, SubEpiPredictWarning

T = np.arange(0, 30, dtype=float)


def context(ydata, n=1, flag1=3, method1=0, thr=0.0, I0=None, onset_fixed=0):
    return sp.Struct(timevect=T, ydata=np.asarray(ydata, dtype=float), I0=float(ydata[0]) if I0 is None else I0,
                     npatches=n, onset_fixed=onset_fixed, onset_thr=thr, flag1=flag1, method1=method1,
                     odeOptions=sp.subepidemicOdeOptions())  # fmt: skip


def logistic_curve(r=0.45, K=4000.0, I0=4.0):
    return sp.simulateSubepidemic(T, I0, 1, 0, 0, 3, [r], [1], [1], [K])[2]


# ---- one candidate ---------------------------------------------------------------


def test_recovers_the_parameters_of_a_noise_free_logistic_curve():
    y = logistic_curve()
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y, 3, 0)
    P, fval, info = sp.fitSubepidemicModel(context(y), z0, LB, UB, isLog, None, 6)
    assert info.success and info.exitflag > 0 and info.nActive == 1
    assert P[0] == pytest.approx(0.45, rel=1e-4) and P[3] == pytest.approx(4000.0, rel=1e-4)
    assert P[1] == 1 and P[2] == 1 and P[4] == 0 and P[5] == 0  # fixed parameters keep their values
    assert fval < 1e-3 * np.sum(y**2) * 1e-6
    assert info.nFree == 2 and info.freeIndex.tolist() == [0, 3] and info.isLog.tolist() == [True, False, False, True, False, False]
    assert info.nLocalRuns == 10 and 1 <= info.nConvergedRuns <= 10 and 0 < info.startAgreement <= 1
    assert info.attempts == 1 and info.atLowerBound.size == 0 and info.atUpperBound.size == 0
    assert info.fallback is None and info.gradient.size == 2 and math.isfinite(info.gradientInfNorm)
    assert info.optimizerOptions.FiniteDifferenceStepSize == 1e-6 and info.useParallel is False
    assert not info.polish.run or info.polish.objectiveAfter <= info.polish.objectiveBefore


@pytest.mark.parametrize("method1, alpha", [(1, None), (3, 0.8), (4, 0.02), (6, None)])
def test_other_estimation_methods_recover_the_curve(method1, alpha):
    mean = logistic_curve()
    rng = np.random.default_rng(12)
    if method1 == 1:
        y = rng.poisson(mean).astype(float)
    elif method1 == 6:
        y = mean + rng.laplace(0, 3, mean.size)
        y[0] = mean[0]
    else:
        dist1 = method1
        y = sp.AddPoissonError(np.cumsum(mean), 1, dist1, alpha, rng=rng)[:, 0]
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y, 3, method1)
    P, fval, info = sp.fitSubepidemicModel(context(y, method1=method1), z0, LB, UB, isLog, None, 8)
    assert info.success
    assert P[0] == pytest.approx(0.45, rel=0.05) and P[3] == pytest.approx(4000.0, rel=0.08)
    value, detail = sp.subepidemicObjective(P, context(y, method1=method1))
    assert value == pytest.approx(fval) and detail.success
    truth = np.array([0.45, 1, 1, 4000.0, P[4], P[5]])
    assert fval <= sp.subepidemicObjective(truth, context(y, method1=method1))[0] + 1e-6  # at least as good as the truth
    if alpha is not None:
        assert 0 < P[4] < 10 * alpha + 1 and info.nFree == 3


def test_nothing_to_estimate_evaluates_the_fixed_vector():
    y = logistic_curve()
    z = np.array([0.45, 1, 1, 4000.0, 0, 0])
    P, fval, info = sp.fitSubepidemicModel(context(y), z, z, z, np.zeros(6, dtype=bool))
    assert info.success and info.nFree == 0 and np.array_equal(P, z) and fval == pytest.approx(0, abs=1e-6)


def test_a_fit_with_no_valid_point_is_reported_not_replaced():
    y = logistic_curve()
    y0 = y.copy()
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y0, 3, 0)
    # initial value 0: every trajectory is identically zero
    P, fval, info = sp.fitSubepidemicModel(context(y0, I0=0.0), z0, LB, UB, isLog, None, 3)
    assert not info.success and np.all(np.isnan(P)) and fval == math.inf
    assert info.fallback is None and info.attempts == 3 and "identically zero" in info.message
    assert info.startAgreement == 1.0  # the penalty plateau "converges" from every start


def test_supplied_start_points_are_used_and_checked():
    y = logistic_curve()
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y, 3, 0)
    good = np.array([[0.4, 1, 1, 5000.0, 0, 0], [0.5, 1, 1, 3000.0, 0, 0]])
    P, fval, info = sp.fitSubepidemicModel(context(y), z0, LB, UB, isLog, good, 0)
    assert info.success and info.nLocalRuns == 2 and P[0] == pytest.approx(0.45, rel=1e-4)
    for args, identifier in (
        ((context(y), z0, LB, UB, isLog, np.zeros((2, 4)), 0), "SubEpiPredict:ParameterDimensionMismatch"),
        ((context(y), z0[:4], LB, UB, isLog), "SubEpiPredict:ParameterDimensionMismatch"),
        ((context(y, method1=2), z0, LB, UB, isLog), "SubEpiPredict:InvalidMethod"),
        (({"timevect": T}, z0, LB, UB, isLog), "SubEpiPredict:FitContext"),
        ((context(y), z0, LB, UB, isLog, None, 2, dict(optimOptions={"Colour": 1})), "SubEpiPredict:FitSettings"),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.fitSubepidemicModel(*args)
        assert identifier_of(e) == identifier


def test_same_seed_same_fit_and_parallel_equals_serial():
    y = np.round(logistic_curve() * (1 + 0.1 * np.cos(T)))
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y, 1, 3)
    c = context(y, flag1=1, method1=3)
    sp.rng(8)
    P1, f1, info1 = sp.fitSubepidemicModel(c, z0, LB, UB, isLog, None, 5)
    sp.rng(8)
    P2, f2, _ = sp.fitSubepidemicModel(c, z0, LB, UB, isLog, None, 5)
    assert np.array_equal(P1, P2) and f1 == f2
    sp.rng(8)
    try:
        P3, f3, info3 = sp.fitSubepidemicModel(c, z0, LB, UB, isLog, None, 5, dict(useParallel=True))
    finally:
        sp.shutdown_pool()
    assert np.array_equal(P1, P3) and f1 == f3 and info1.useParallel is False
    assert info3.useParallel in (True, False)  # False when no pool could be started


def test_estimates_on_a_bound_are_listed():
    y = logistic_curve(K=4000.0)
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(1, y, 3, 0, dict(KUB=1500.0))
    P, fval, info = sp.fitSubepidemicModel(context(y), z0, LB, UB, isLog, None, 4, dict(KUB=1500.0))
    assert info.success and P[3] == pytest.approx(1500.0) and info.atUpperBound.tolist() == [3]


# ---- threshold grid ------------------------------------------------------------------


def test_threshold_grid(covid_series):
    grid = sp.onsetThresholdGrid(covid_series, 7)
    assert grid.size == 75 and grid[0] == 0 and grid[1] == pytest.approx(0.5)
    cum = np.cumsum(sp.movmean(covid_series, 7))
    assert np.allclose(np.diff(grid[1:]), (cum[-1] - cum[0]) / 74)
    # the thresholds of the four ranked models in the README figure
    for value in (1140.5536, 27361.7865, 28501.8401):
        assert np.min(np.abs(grid - value)) < 1e-4


# ---- screen and bootstrap on a synthetic two-sub-epidemic series -----------------------

TRUTH = dict(r=[0.5, 0.6], K=[3000.0, 2000.0], thr=1500.0, I0=5.0)


@pytest.fixture(scope="module")
def synthetic():
    curve = sp.simulateSubepidemic(T, TRUTH["I0"], 2, 0, TRUTH["thr"], 3, TRUTH["r"], [1, 1], [1, 1], TRUTH["K"])[2]
    return np.column_stack([T, curve])


@pytest.fixture(scope="module")
def screen(synthetic, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("screen"))
    sp.rng(1)
    result = sp.fittingModifiedLogisticFunctionPatchABC(
        "daily-test-cases-X-state-1-01-01-2020.txt", synthetic, 1, np.arange(1, 31), 5, 3, 4,
        method1=0, dist1=0, npatches_fixed=2, onset_fixed=0, smoothfactor1=1, calibrationperiod1=30,
        output_dir=out, plot=False, verbose=False,
    )  # fmt: skip
    return result, out


@pytest.mark.slow
def test_screen_finds_the_two_sub_epidemics(screen, synthetic):
    (RMSES, PS, npatches, onset_thr, P, info), out = screen
    grid = sp.onsetThresholdGrid(synthetic[:, 1], 1)
    nearest = grid[np.argmin(np.abs(grid - TRUTH["thr"]))]
    assert npatches == 2 and onset_thr == pytest.approx(nearest)
    assert np.allclose(P[0:2], TRUTH["r"], rtol=0.02) and np.allclose(P[6:8], TRUTH["K"], rtol=0.02)
    # ranking: sorted by AICc, one row per retained candidate, parameters stored row by row
    assert np.all(np.diff(RMSES[:, 2]) >= 0) and PS.shape == (RMSES.shape[0], 10)
    assert np.array_equal(PS[0], P)
    one = RMSES[RMSES[:, 0] == 1]
    assert one.shape[0] == 1 and one[0, 1] == 0 and np.all(PS[RMSES[:, 0] == 1, 6:] == 0)  # unused columns are zero
    assert np.allclose(info.relativelik_i, np.exp((RMSES[0, 2] - RMSES[:, 2]) / 2)) and info.relativelik_i[0] == 1
    # every candidate kept has all its sub-epidemics active; thresholds are not repeated
    two = RMSES[RMSES[:, 0] == 2, 1]
    assert np.unique(two).size == two.size and 5.0 not in two.tolist()  # 5 = I0 is met at time zero, like 0
    d = info.thresholdDiagnostics
    assert d.shape[1] == 5 and np.sum(d[:, 3] == 1) == RMSES.shape[0]
    assert np.all(d[d[:, 3] == 1, 2] == d[d[:, 3] == 1, 0])
    assert info.fitDiagnostics.shape == (d.shape[0], 11) and len(info.candidateFitInfo) == d.shape[0]
    assert np.allclose(info.bestfit, synthetic[:, 1], rtol=0.05, atol=2)
    # extra start points for the first fit of each sweep: 3 times numstartpoints, plus the seed
    assert info.fitDiagnostics[0, 8] == 13 and info.fitDiagnostics[1, 8] == 13
    assert set(info.fitDiagnostics[2:, 8].tolist()) <= {5.0, 6.0, 13.0}


@pytest.mark.slow
def test_screen_file_can_be_read_back(screen):
    (RMSES, PS, npatches, onset_thr, P, info), out = screen
    name = ("ABC-ensem-npatchesfixed-2-onsetfixed-0-smoothing-1-daily-test-cases-X-state-1-01-01-2020"
            "-flag1-3-method-0-dist-0-calibrationperiod-30.mat")  # fmt: skip
    assert os.path.basename(info.file) == name and os.path.isfile(info.file)
    saved = sp.load_mat(info.file)
    assert np.allclose(saved.RMSES, RMSES) and np.allclose(saved.PS, PS)
    assert saved.data1.shape == (30, 2) and float(saved.I0[0, 0]) == 5.0


@pytest.mark.slow
def test_bootstrap_of_the_best_model(screen, synthetic):
    (RMSES, PS, _, _, _, _), _ = screen
    sp.rng(2)
    out = sp.fittingModifiedLogisticFunctionPatchMultiple(
        RMSES, PS, synthetic, 1, np.arange(1, 31), 12, 3, 4, 1,
        method1=0, dist1=0, onset_fixed=0, smoothfactor1=1, plot=False, verbose=False,
    )  # fmt: skip
    Phatss, npatches, onset_thr, curves, bestfit, data1, P0, AICc_best, factor1, d, diag = out
    assert Phatss.shape == (12, 10) and curves.shape == (30, 12) and npatches == 2
    assert np.array_equal(P0, PS[0]) and AICc_best == RMSES[0, 2] and onset_thr == RMSES[0, 1] and d == 0
    assert np.all(curves[0] == synthetic[0, 1])  # the first observation is not perturbed
    # normal errors: the noise scale is the residual standard deviation of the fit
    resid = bestfit - synthetic[:, 1]
    assert factor1 == pytest.approx(np.sqrt(np.sum(resid**2) / (30 - 6)))
    assert np.allclose(np.median(Phatss[:, [0, 1, 6, 7]], axis=0), [0.5, 0.6, 3000, 2000], rtol=0.03)
    assert diag.success.sum() + diag.usedFallback.sum() == 12 and np.all(diag.nActive == 2)
    assert diag.searchLB.size == 10 and diag.identifiability.M == 12
    assert diag.identifiability.name[:2] == ["r_1", "r_2"] and diag.notIdentified == []
    # rank 2 can be bootstrapped too; an invalid ranked vector is an error
    bad = PS.copy()
    bad[0, 0] = np.nan
    with pytest.raises(SubEpiPredictError) as e:
        sp.fittingModifiedLogisticFunctionPatchMultiple(
            RMSES, bad, synthetic, 1, np.arange(1, 31), 2, 3, 4, 1,
            method1=0, dist1=0, onset_fixed=0, smoothfactor1=1, plot=False, verbose=False,
        )  # fmt: skip
    assert identifier_of(e) == "SubEpiPredict:RankedFitInvalid"


@pytest.mark.slow
def test_models_too_large_for_the_series_are_left_out():
    y = logistic_curve()[:9]
    data1 = np.column_stack([T[:9], y])
    with pytest.warns(SubEpiPredictWarning, match="CandidateTooLarge"):
        RMSES, PS, npatches, _, _, info = sp.fittingModifiedLogisticFunctionPatchABC(
            "x.txt", data1, 1, np.arange(1, 10), 5, 1, 1, method1=0, dist1=0, npatches_fixed=3, onset_fixed=1,
            smoothfactor1=1, calibrationperiod1=9, save=False, plot=False, verbose=False,
        )  # fmt: skip
    # 3 GLM sub-epidemics with normal errors have 10 parameters: too many for 9 points
    assert info.candidatesTooLarge.tolist() == [[3, 10, 9]] and set(RMSES[:, 0].tolist()) <= {1.0, 2.0}


# ---- the README figure of the MATLAB toolbox -------------------------------------------


@pytest.mark.slow
def test_top_ranked_model_of_the_quick_start_matches_the_matlab_figures(covid_series):
    """docs/images/parameters.png and ranking_profiles.png show ranked model 1 of the
    quick start: C_thr = 1140.5536, r = 0.561 and 4.63, p = 0.826 and 0.659,
    K = 80400 and 24800, AICc = -1.09196e+06."""
    y = covid_series
    thr = sp.onsetThresholdGrid(y, 7)[2]
    assert thr == pytest.approx(1140.5536, abs=1e-4)
    LB, UB, isLog, z0 = sp.subepidemicSearchSpace(2, y, 1, 3)
    c = sp.Struct(timevect=np.arange(75.0), ydata=y, I0=y[0], npatches=2, onset_fixed=0, onset_thr=thr, flag1=1,
                  method1=3, odeOptions=sp.subepidemicOdeOptions())  # fmt: skip
    sp.rng(1)
    m = sp.subepidemicUnitMap("make", dict(z0=z0, LB=LB, UB=UB, isLog=isLog))
    starts = sp.subepidemicStartPoints(m, 60, 0, z0)
    P, fval, info = sp.fitSubepidemicModel(c, z0, LB, UB, isLog, starts, 60)
    assert info.success and info.nActive == 2
    AICc = sp.getAICc(3, 3, 2, 1, 1, fval, 75, 0)[0]
    assert "%.6g" % AICc == "-1.09196e+06" and AICc == pytest.approx(-1091961.08, abs=0.05)
    assert P[0] == pytest.approx(0.561, abs=0.002) and P[1] == pytest.approx(4.63, abs=0.02)
    assert P[2] == pytest.approx(0.826, abs=0.002) and P[3] == pytest.approx(0.659, abs=0.002)
    assert P[6] == pytest.approx(80400, abs=200) and P[7] == pytest.approx(24800, abs=100)
