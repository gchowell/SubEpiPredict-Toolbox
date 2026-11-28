"""Parameter identifiability: PII, bound hits, Monte Carlo error, reports."""

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError


def test_sci_interval():
    assert sp.sci_interval(10.0, 1000.0) == pytest.approx(2.0)
    assert sp.sci_interval(5.0, 5.0) == 0
    assert sp.sci_interval(0.0, 3.0) == pytest.approx(6.0)  # the cap: not bounded away from zero
    assert np.isnan(sp.sci_interval(np.nan, 3.0))
    out = sp.sci_interval(np.array([1.0, 0.0, 2.0]), np.array([10.0, 0.0, 2.0]))
    assert out.tolist() == [1.0, 0.0, 0.0]
    assert sp.sci_interval(1e-9, 1.0) == pytest.approx(6.0)  # floored at 1e-6 of the upper limit


def test_pii_definitions():
    lo = np.array([1.0, 0.2, 0.0, 5.0])
    hi = np.array([100.0, 0.4, 3.0, 5.0])
    LB = np.array([0.1, 0.0, 0.0, 5.0])
    UB = np.array([1000.0, 1.0, 10.0, 5.0])
    PII, rangeFrac, scale, PIIbox, ratio = sp.subepidemicPII(lo, hi, LB, UB)
    assert PII[0] == pytest.approx(2) and PII[1] == pytest.approx(np.log10(2)) and np.isnan(PII[2]) and PII[3] == 0
    assert np.allclose(rangeFrac[:3], [99 / 999.9, 0.2, 0.3]) and np.isnan(rangeFrac[3])
    assert scale.tolist() == [1, 0, 0, 1]
    assert PIIbox[0] == pytest.approx(4) and np.isnan(PIIbox[1]) and PIIbox[3] == 0
    assert ratio[0] == pytest.approx(0.5) and ratio[1] == pytest.approx(0.2) and ratio[2] == pytest.approx(0.3)
    assert np.isnan(ratio[3])  # fixed parameter
    # the user's scale replaces the automatic one where given
    _, _, scale2, _, ratio2 = sp.subepidemicPII(lo, hi, LB, UB, [0, np.nan])
    assert scale2.tolist() == [0, 0, 0, 1] and ratio2[0] == pytest.approx(99 / 999.9)
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicPII(lo, hi, LB, UB, [2])
    assert identifier_of(e) == "SubEpiPredict:PIIScale"
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicPII(lo, hi[:2], LB, UB)
    assert identifier_of(e) == "SubEpiPredict:PIISize"


def test_bound_hit_fractions():
    P = np.array([[1.0, 5.0, 0.0], [1.0005, 9.995, 0.0], [3.0, 10.0, 0.0], [np.nan, 7.0, 0.0]])
    LB = [1.0, 5.0, 0.0]
    UB = [10.0, 10.0, 0.0]
    fracLB, fracUB, hitLB, hitUB = sp.subepidemicBoundHitFractions(P, LB, UB)
    assert fracLB[0] == pytest.approx(2 / 3) and fracUB[0] == 0  # NaN is not counted
    assert fracLB[1] == 0.25 and fracUB[1] == 0.5  # 9.995 is within 0.1% of 10
    assert np.isnan(fracLB[2]) and np.isnan(fracUB[2])  # fixed column
    assert hitLB[:, 0].tolist() == [True, True, False, False]
    assert sp.subepidemicBoundHitFractions(P, LB, UB, 1e-6)[0][0] == pytest.approx(1 / 3)
    # a flipped box is put right; a wrong size is an error
    assert sp.subepidemicBoundHitFractions(P, UB, LB)[0][1] == 0.25
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicBoundHitFractions(P, [1.0], [2.0])
    assert identifier_of(e) == "SubEpiPredict:BoundHitSizeMismatch"


def test_quantile_mcse_is_repeatable_and_does_not_touch_the_stream():
    rng = np.random.default_rng(4)
    P = np.column_stack([rng.lognormal(0, 0.5, 60), np.full(60, 3.0), rng.normal(0, 1, 60)])
    sp.rng(9)
    before = sp.get_rng().bit_generator.state["state"]["pos"]
    a = sp.subepidemicQuantileMCSE(P)
    assert sp.get_rng().bit_generator.state["state"]["pos"] == before
    b = sp.subepidemicQuantileMCSE(P)
    for u, v in zip(a, b):
        assert np.array_equal(u, v, equal_nan=True)
    mcseLB, mcseUB, mcsePII = a
    assert mcseLB[0] > 0 and mcseUB[0] > 0 and mcsePII[0] > 0
    assert np.isnan(mcseLB[1]) and np.isnan(mcsePII[1])  # no spread
    assert mcseLB[2] > 0 and np.isnan(mcsePII[2])  # lower quantile below zero: no log ratio
    assert all(np.all(np.isnan(v)) for v in sp.subepidemicQuantileMCSE(P[:1]))


def bootstrap_matrix(M=50):
    rng = np.random.default_rng(1)
    r = rng.lognormal(np.log(0.4), 0.1, (M, 2))
    p = np.clip(rng.normal(0.8, 0.05, (M, 2)), 0, 1)
    a = np.ones((M, 2))
    K = rng.lognormal(np.log(5e4), 0.3, (M, 2))
    K[:10, 1] = 1e10  # a fifth of the refits at the upper limit of K_2
    alpha = rng.lognormal(0, 0.2, (M, 1))
    return np.hstack([r, p, a, K, alpha, np.zeros((M, 1))])


LB = np.array([0.005, 0.005, 0, 0, 1, 1, 20, 20, 1e-8, 0])
UB = np.array([500, 500, 1, 1, 1, 1, 1e10, 1e10, 1e5, 0])


def test_identifiability_report():
    P = bootstrap_matrix()
    P0 = np.median(P, axis=0)
    pii = sp.subepidemicIdentifiability(P, LB, UB, 2, P0)
    assert pii.name == ["r_1", "r_2", "p_1", "p_2", "a_1", "a_2", "K_1", "K_2", "alpha"]  # d is fixed: not listed
    assert pii.type[-1] == "alpha" and np.isnan(pii.subepidemic[-1]) and pii.subepidemic[1] == 2
    assert pii.column.tolist() == [0, 1, 2, 3, 4, 5, 6, 7, 8] and pii.M == 50 and pii.B2 == 500
    assert np.allclose(pii.CI_LB, sp.quantile(P, 0.025, axis=0)[:9]) and np.allclose(pii.estimate, P0[:9])
    assert pii.fixed.tolist() == [False, False, False, False, True, True, False, False, False]
    assert pii.PII_scale.tolist() == [1, 1, 0, 0, 1, 1, 1, 1, 1]
    assert pii.fracUB[7] == pytest.approx(0.2) and pii.fracUB[6] == 0
    assert pii.PII_ratio[7] > 0.5 > pii.PII_ratio[6]  # K_2 reaches its limit: interval set by the bounds
    assert pii.PII_box[6] == pytest.approx(np.log10(1e10 / 20))
    # user scale for r: range fraction
    pii2 = sp.subepidemicIdentifiability(P, LB, UB, 2, P0, [0, np.nan, np.nan, np.nan, np.nan, np.nan])
    assert pii2.PII_scale[:2].tolist() == [0, 0] and pii2.PII_ratio[0] == pytest.approx(pii.PII_rangeFrac[0])
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicIdentifiability(P[:, :8], LB, UB, 2)


def test_bound_summary():
    P = bootstrap_matrix()
    isLog = np.array([1, 1, 0, 0, 0, 0, 1, 1, 1, 0], dtype=bool)
    report, lines, nRows = sp.subepidemicBoundSummary(P, LB, UB, isLog, 2)
    K2 = [item for item in report if item.name == "K2"]
    assert len(K2) == 1 and K2[0].side == "upper" and K2[0].count == 10 and K2[0].total == 50
    assert K2[0].fraction == 0.2 and K2[0].setting == "KUB" and K2[0].isSearchLimit and K2[0].index == 7
    assert "K2 at its upper limit (1e+10) in 10 of 50" in lines
    assert nRows >= 10
    # p at its upper bound is reported but is not a search limit
    p_items = [item for item in report if item.name.startswith("p")]
    assert all(not item.isSearchLimit and item.setting in ("pLB", "pUB") for item in p_items)
    # a single vector, and rows with NaN are ignored
    one, _, n1 = sp.subepidemicBoundSummary(UB, LB, UB, isLog, 2)
    assert {item.name for item in one} >= {"r1", "r2", "K1", "K2", "alpha"} and n1 == 1
    assert sp.subepidemicBoundSummary(np.full((3, 10), np.nan), LB, UB, isLog, 2) == ([], [], 0)


def test_print_and_write_pii(tmp_path, capsys):
    P = bootstrap_matrix()
    pii = sp.subepidemicIdentifiability(P, LB, UB, 2, np.median(P, axis=0))
    sp.subepidemicPrintPII(pii, "PII of the test")
    out = capsys.readouterr().out
    assert "PII of the test (50 bootstrap refits):" in out and "fixed" in out and "K_2" in out
    one = tmp_path / "one.csv"
    sp.subepidemicWritePII(str(one), pii)
    lines = one.read_text().splitlines()
    assert lines[0].startswith("parameter,subepidemic,estimate,median,CI_LB,CI_UB,PII,") and lines[0].endswith(",fixed")
    assert len(lines) == 10 and lines[1].startswith("r_1,1,") and lines[-1].startswith("alpha,NaN,")
    assert lines[5].endswith(",1")  # a_1 is fixed
    both = tmp_path / "both.csv"
    sp.subepidemicWritePII(str(both), [pii, pii], [1, 2])
    lines = both.read_text().splitlines()
    assert lines[0].startswith("ranked_model,n_subepidemics,parameter") and len(lines) == 19
    assert lines[1].startswith("1,2,r_1,1,") and lines[10].startswith("2,2,r_1,1,")
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicWritePII(str(both), [pii, pii], [1])
    assert identifier_of(e) == "SubEpiPredict:PIIRanks"
