"""Quantiles, weighted interval score, performance metrics, baseline, ensemble weights."""

import os

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError, SubEpiPredictWarning
from subepipredict.scores import WIS_ALPHAS

RNG = np.random.default_rng(3)
CURVES = RNG.poisson(80, (12, 300)).astype(float)
P = [0.0, 0.01, 0.025, 0.25, 0.5, 0.9, 0.975, 1.0]

# ---- quantiles ------------------------------------------------------------------


def test_matlab_quantile_definition():
    # quantile([1 2 3 4], p) in MATLAB
    assert sp.quantile([1, 2, 3, 4], 0.5) == 2.5
    assert np.allclose(sp.quantile([1, 2, 3, 4], [0.1, 0.25, 0.3, 0.75, 0.95]), [1, 1.5, 1.7, 3.5, 4])
    assert np.allclose(sp.quantile([[1, 10], [3, 30], [2, 20]], 0.5, axis=0), [2, 20])
    assert sp.quantile([1, np.nan, 3], 0.5) == 2


def test_weighted_quantile_with_equal_weights_is_the_plain_quantile():
    q = sp.subepidemicWeightedQuantile(CURVES, P)
    assert q.shape == (12, len(P))
    assert np.allclose(q, sp.quantile(CURVES, P, axis=1).T, atol=1e-10)
    assert np.allclose(sp.subepidemicWeightedQuantile(CURVES, P, np.full(300, 7.0)), q)


def test_weighted_quantile_definition_and_invariances():
    assert sp.subepidemicWeightedQuantile([[1, 2, 3, 4]], [0.5], [1, 1, 1, 5])[0, 0] == pytest.approx(3.5)
    w = RNG.uniform(0.1, 2, 300)
    q = sp.subepidemicWeightedQuantile(CURVES, P, w)
    order = RNG.permutation(300)
    assert np.allclose(sp.subepidemicWeightedQuantile(CURVES[:, order], P, w[order]), q)  # order of the curves
    assert np.allclose(sp.subepidemicWeightedQuantile(CURVES, P, 10 * w), q)  # scale of the weights
    # a curve of weight 2 counts as two copies of weight 1
    x = np.array([[5.0, 1.0, 9.0, 4.0]])
    doubled = sp.subepidemicWeightedQuantile(np.array([[5.0, 1.0, 9.0, 4.0, 9.0]]), P)
    assert np.allclose(sp.subepidemicWeightedQuantile(x, P, [1, 1, 2, 1]), doubled)
    # zero weights and NaN values are left out; nothing left gives NaN
    q0 = sp.subepidemicWeightedQuantile([[1.0, 50.0, 3.0, np.nan]], [0.5], [1, 0, 1, 1])
    assert q0[0, 0] == 2.0
    assert np.isnan(sp.subepidemicWeightedQuantile([[np.nan, 2.0]], [0.5], [1, 0])[0, 0])


def test_weighted_quantile_errors():
    for args in (([[1, 2]], [1.5]), ([[1, 2]], [0.5], [1]), ([[1, 2]], [0.5], [1, -1]), ([[1, 2]], [])):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicWeightedQuantile(*args)
        assert identifier_of(e) == "SubEpiPredict:WeightedQuantile"


def test_compute_quantiles_splits_calibration_and_forecast():
    data1 = np.column_stack([np.arange(9.0), np.ones(9)])
    curves = CURVES - 75  # some values below zero
    qc, qf, alphas = sp.computeQuantiles(data1, curves, 3)
    assert qc.shape == (9, 23) and qf.shape == (3, 23) and alphas.size == 23
    full = sp.quantile(curves, alphas, axis=1).T
    assert np.allclose(qc, np.where(full[:9] >= 0, full[:9], 0)) and np.allclose(qf, np.where(full[9:] >= 0, full[9:], 0))
    assert np.all(qc >= 0)
    qc0, qf0, _ = sp.computeQuantiles(data1, curves, 0)
    assert qf0.shape == (0, 23) and np.array_equal(qc0, qc)


# ---- weighted interval score ------------------------------------------------------


def manual_wis(y, samples):
    total = 0.5 * abs(y - np.median(samples))
    for alpha in WIS_ALPHAS:
        L, U = sp.quantile(samples, [alpha / 2, 1 - alpha / 2])
        IS = (U - L) + (2 / alpha) * (L - y) * (y < L) + (2 / alpha) * (y - U) * (y > U)
        total += (alpha / 2) * IS
    return total / (len(WIS_ALPHAS) + 0.5)


def test_wis_matches_its_definition():
    obs = np.array([70.0, 82, 95, 60, 81, 120, 79, 80, 77, 90, 150, 30])
    data1 = np.column_stack([np.arange(9.0), obs[:9]])
    latest = np.column_stack([np.arange(12.0), obs])
    WISC, WISFS, qc, qf = sp.computeWIS(data1, latest, CURVES, 3)
    each = np.array([manual_wis(obs[i], CURVES[i]) for i in range(12)])
    assert WISC == pytest.approx(each[:9].mean(), rel=1e-12)
    assert WISFS[:, 0].tolist() == [1, 2, 3]
    assert np.allclose(WISFS[:, 1], [each[9], each[9:11].mean(), each[9:12].mean()])
    assert qc.shape == (9, 23) and qf.shape == (3, 23)
    # no forecast period
    WISC0, WISFS0, _, _ = sp.computeWIS(data1, data1, CURVES[:9], 0)
    assert WISC0 == pytest.approx(WISC) and WISFS0.shape == (0, 2)
    # a point forecast: WIS reduces to a multiple of the absolute error
    point = np.full((9, 5), 100.0)
    w, _, _, _ = sp.computeWIS(data1, data1, point, 0)
    expected = np.mean((0.5 * np.abs(obs[:9] - 100) + np.sum([np.abs(obs[:9] - 100) for _ in WIS_ALPHAS], axis=0)) / 11.5)
    assert w == pytest.approx(expected)


def test_performance_metrics():
    obs = np.array([70.0, 82, 95, 60, 81, 120, 79, 80, 77, 90, 150, 30])
    data1 = np.column_stack([np.arange(9.0), obs[:9]])
    latest = np.column_stack([np.arange(12.0), obs])
    traj = CURVES[:, :40] + 0.5
    r = sp.computeforecastperformance(data1, latest, traj, CURVES, 3)
    RMSEc, MSEc, MAEc, PIc, MISc, RMSEF, MSEF, MAEF, PIF, MISF = r
    med = np.median(traj, axis=1)
    L = sp.quantile(CURVES, 0.025, axis=1)
    U = sp.quantile(CURVES, 0.975, axis=1)
    err = obs - med
    assert MSEc == pytest.approx(np.mean(err[:9] ** 2)) and RMSEc == pytest.approx(np.sqrt(MSEc))
    assert MAEc == pytest.approx(np.mean(np.abs(err[:9])))
    inside = (obs >= L) & (obs <= U)
    assert PIc == pytest.approx(100 * inside[:9].mean())
    IS = (U - L) + 40 * (L - obs) * (obs < L) + 40 * (obs - U) * (obs > U)
    assert MISc == pytest.approx(IS[:9].mean())
    for h in (1, 2, 3):
        assert MAEF[h - 1].tolist() == pytest.approx([h, np.mean(np.abs(err[9 : 9 + h]))])
        assert MSEF[h - 1, 1] == pytest.approx(np.mean(err[9 : 9 + h] ** 2))
        assert RMSEF[h - 1, 1] == pytest.approx(np.sqrt(MSEF[h - 1, 1]))
        assert PIF[h - 1, 1] == pytest.approx(100 * inside[9 : 9 + h].mean())
        assert MISF[h - 1, 1] == pytest.approx(IS[9 : 9 + h].mean())
    r0 = sp.computeforecastperformance(data1, data1, traj[:9], CURVES[:9], 0)
    assert r0[2] == pytest.approx(MAEc) and r0[7].shape == (0, 2)


def test_weights_enter_the_scores():
    obs = np.full(12, 80.0)
    data1 = np.column_stack([np.arange(12.0), obs])
    w = np.zeros(300)
    w[:150] = 1  # only the first half of the curves counts
    a = sp.computeWIS(data1, data1, CURVES, 0, w)[0]
    b = sp.computeWIS(data1, data1, CURVES[:, :150], 0)[0]
    assert a == pytest.approx(b, rel=1e-12)
    pa = sp.computeforecastperformance(data1, data1, CURVES, CURVES, 0, w, w)
    pb = sp.computeforecastperformance(data1, data1, CURVES[:, :150], CURVES[:, :150], 0)
    assert np.allclose(pa[:5], pb[:5])


# ---- baseline (the numbers of the README) ----------------------------------------


def test_baseline_of_the_quick_start(covid_series, input_dir):
    latest = np.loadtxt(os.path.join(input_dir, "cumulative-daily-coronavirus-deaths-USA-05-09-2022.txt"), ndmin=2)[:, 51]
    later = np.diff(latest)[74:78]
    assert later.tolist() == [1589, 1735, 1773, 1665]  # the four targets after 05-11-2020
    B = sp.subepidemicBaselineForecast(covid_series, 4, 7)
    assert B.level == pytest.approx(1675.7142857142858) and round(B.level, 2) == 1675.71
    assert np.all(B.point == B.level) and B.window == 7
    assert np.mean(np.abs(later - B.point)) == pytest.approx(63.5, abs=1e-8)
    assert B.nErrors.tolist() == [68, 67, 66, 65] and not B.borrowed.any()
    assert B.samples.shape == (4, 200) and np.all(B.samples >= 0)
    assert np.allclose(np.median(B.samples, axis=1), B.level)  # symmetric errors keep the median at the level


def test_baseline_rules():
    y = np.array([3.0, 5, 4, 8, 7, 9, 12, 10])
    B = sp.subepidemicBaselineForecast(y, 2, 1)
    assert B.level == 10 and B.nErrors.tolist() == [7, 6]
    e1 = np.diff(y)
    both = np.sort(np.concatenate([e1, -e1]))
    assert np.allclose(B.samples[0], np.maximum(10 + sp.quantile(both, (np.arange(200) + 0.5) / 200), 0))
    # too few errors at a long horizon: those of the longest shorter horizon are borrowed
    B = sp.subepidemicBaselineForecast(y, 5, 1)
    assert B.nErrors.tolist() == [7, 6, 5, 5, 5] and B.borrowed.tolist() == [False, False, False, True, True]
    assert np.array_equal(B.samples[3], B.samples[2])
    # no horizon with five errors: the uncertainty is not defined
    B = sp.subepidemicBaselineForecast(y[:4], 2, 1)
    assert np.all(np.isnan(B.samples)) and B.level == 8
    assert sp.subepidemicBaselineForecast(y, 1, 50).window == 8
    for args in (([], 1, 1), (y, 0, 1), (y, 1, 0), (y, 1.5, 1)):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicBaselineForecast(*args)
        assert identifier_of(e) == "SubEpiPredict:Baseline"


# ---- ensemble weights and evidence ratios ------------------------------------------


def test_ensemble_weights():
    AICc = np.array([-1091961.08, -1091959.62, -1091957.92])
    assert np.allclose(sp.subepidemicEnsembleWeights(-1, AICc), [1 / 3] * 3)
    w = sp.subepidemicEnsembleWeights(1, AICc)
    rel = np.exp(-(AICc - AICc.min()) / 2)
    assert np.allclose(w, rel / rel.sum()) and w[0] > w[1] > w[2] and w.sum() == pytest.approx(1)
    with pytest.warns(SubEpiPredictWarning, match="WeightCodeZero"):
        assert np.allclose(sp.subepidemicEnsembleWeights(0, AICc), w)
    assert np.allclose(sp.subepidemicEnsembleWeights(2, AICc, [1.0, 2.0, 4.0]), np.array([4, 2, 1]) / 7)
    assert np.allclose(sp.subepidemicEnsembleWeights(3, AICc, None, [1.0, 1.0, 2.0]), np.array([2, 2, 1]) / 5)
    for args in ((5, AICc), (1, [1.0, np.inf]), (2, AICc, [1.0, 2.0]), (2, AICc, [1.0, 0.0, 2.0]), (3, AICc)):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicEnsembleWeights(*args)
        assert identifier_of(e) == "SubEpiPredict:EnsembleWeights"


def test_evidence_ratios():
    AICc = np.array([100.0, 101.4, 110.0, 130.0, 140.0])
    ratio, weak = sp.subepidemicEvidenceRatios(AICc, 2, 1)
    assert np.allclose(ratio, [1, np.exp(0.7)]) and weak.size == 0
    with pytest.warns(SubEpiPredictWarning, match="WeakTopModel") as record:
        ratio, weak = sp.subepidemicEvidenceRatios(AICc, 4, 52)
    assert weak.tolist() == [3, 4] and ratio[2] == pytest.approx(np.exp(5))
    assert "ranked models 3, 4 have evidence ratios of 148, 3.3e+06" in str(record[0].message)
    ratio, weak = sp.subepidemicEvidenceRatios(AICc, 4, 52, np.inf)
    assert weak.size == 0
    ratio, weak = sp.subepidemicEvidenceRatios(AICc, 9, 52, np.inf)
    assert ratio.size == 5
    assert sp.subepidemicEvidenceRatios([], 3, 1)[0].size == 0
