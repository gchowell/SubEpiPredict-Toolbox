"""Reproduction number, doubling times, and the small helpers."""

import math

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError
from subepipredict.utils import factor, mat2str, num2str

T = np.arange(0, 60, dtype=float)

# ---- reproduction number ---------------------------------------------------------


def test_constant_incidence_gives_rt_one_after_the_blanked_start():
    Rs, ts, coverage = sp.get_Rt(T, np.full(60, 100.0), 1, 11.7, 3.7)
    assert ts[0] == 0 and ts[-1] == pytest.approx(59) and ts.size == 5901 and Rs.size == 5900
    assert np.all(np.isnan(Rs[coverage < 0.99 - 1e-9]))
    first = int(np.flatnonzero(np.isfinite(Rs))[0])
    assert 16 < ts[first + 1] < 17.5 and coverage[first] >= 0.99 - 1e-9
    assert np.nanmax(np.abs(Rs[first:] - 1)) < 0.011  # error of at most about 1/coverage - 1
    assert np.all(np.diff(coverage) >= 0) and coverage[-1] == pytest.approx(1)


def test_without_blanking_early_values_are_far_too_large():
    Rs, ts, _ = sp.get_Rt(T, np.full(60, 100.0), 1, 11.7, 3.7, 0)
    assert Rs[499] > 1e5  # day 5: hundreds of thousands
    assert Rs[999] == pytest.approx(5.2, abs=0.1)  # day 10


def test_exponential_growth_against_the_euler_lotka_values():
    r = 0.1
    inc = 100 * np.exp(r * T)
    Rs, ts, _ = sp.get_Rt(T, inc, 3, 5.0, None, 0)  # delta interval: R = exp(r*mean)
    assert np.allclose(Rs[1000:], math.exp(0.5), rtol=1e-3)
    assert np.all(np.isnan(Rs[:499]))  # before one interval has passed there is no infector
    Rs, _, _ = sp.get_Rt(T, inc, 2, 5.0)  # exponential interval: R = 1 + r*mean
    assert Rs[-1] == pytest.approx(1.5, rel=2e-3)
    b = 3.7 / 11.7
    a = 11.7 / b
    Rs, _, _ = sp.get_Rt(T, inc, 1, 11.7, 3.7)  # gamma interval: R = (1 + r*b)^a
    assert Rs[-1] == pytest.approx((1 + r * b) ** a, rel=2e-3)


def test_rt_special_cases_and_errors():
    Rs, ts, cov = sp.get_Rt([3.0], [5.0], 1, 4, 2)
    assert Rs.size == 0 and cov.size == 0 and ts.tolist() == [3.0]
    # a mean that is not a multiple of the grid step
    Rs, _, _ = sp.get_Rt(T, np.full(60, 10.0), 3, 2.337, None, 0)
    assert Rs[-1] == pytest.approx(1)
    for args in (
        (T, np.ones(60), 4, 5),
        (T, np.ones(60), 1, 5),  # the gamma interval needs a variance
        (T, np.ones(60), 1, -5, 2),
        (T, np.ones(59), 1, 5, 2),
        (T, np.ones(60), 1, 5, 2, 1.5),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.get_Rt(*args)
        assert identifier_of(e) == "SubEpiPredict:InvalidRtSettings"


# ---- doubling times -----------------------------------------------------------------


def test_doubling_times_of_a_curve_with_a_constant_doubling_time():
    # cumulative count doubles every 3 time units: C(t) = 2^(t/3)
    t = np.arange(0, 40.0)
    cum = 2 ** (t / 3)
    inc = np.concatenate(([cum[0]], np.diff(cum)))
    tds, C0, curve, doublingtimes = sp.getDoublingTimeCurve(inc, 1)
    assert C0 == 1 and curve.shape == (40, 2) and np.allclose(curve[:, 1], cum)
    # 13 doublings fit in the curve; the one that ends at its last point is not counted
    assert doublingtimes.size == 12 and np.allclose(doublingtimes, 3, atol=1e-9)
    assert np.allclose(tds[1:, 1] / tds[:-1, 1], 2, rtol=1e-9)
    # the time step scales the times
    assert np.allclose(sp.getDoublingTimeCurve(inc, 7)[3], 7 * doublingtimes, atol=1e-9)


def test_doubling_times_special_cases():
    assert sp.getDoublingTimeCurve([5.0], 1)[3].size == 0
    tds, C0, curve, d = sp.getDoublingTimeCurve([10.0, 1, 1, 1], 1)  # never doubles
    assert d.size == 0 and C0 == 10 and curve.shape == (0, 2)
    # cumulative 1, 2, 4, 8, 8, 8: the last level is reached at the end of the plateau and dropped
    tds, C0, curve, d = sp.getDoublingTimeCurve([1.0, 1, 2, 4, 0, 0], 1)
    assert np.allclose(d, [1, 1]) and tds.tolist() == [[0, 1], [1, 2], [2, 4]]


# ---- helpers ---------------------------------------------------------------------------


def test_moving_average_like_smoothdata():
    y = np.array([1.0, 2, 3, 4, 5, 6])
    assert np.allclose(sp.movmean(y, 3), [1.5, 2, 3, 4, 5, 5.5])
    assert np.allclose(sp.movmean(y, 7), [2.5, 3, 3.5, 3.5, 4, 4.5])
    assert np.allclose(sp.movmean(y, 4), [1.5, 2, 2.5, 3.5, 4.5, 5])  # even span: two back, one ahead
    assert np.array_equal(sp.movmean(y, 1), y)


def test_dates():
    assert sp.datenum("05-11-2020") == 737922 == sp.datenum([2020, 5, 11]) == sp.datenum(2020, 5, 11)
    assert sp.datenum(2020, 2, 27) + 74 == sp.datenum("05-11-2020")
    assert sp.datevec(737922) == [2020, 5, 11]
    assert sp.datestr(737922 + 30) == "06-10-2020" and sp.datestr(737922, "mm-dd") == "05-11"
    assert sp.datestr(737922, "mm-dd-yy") == "05-11-20"


def test_text_of_numbers():
    assert num2str(1140.5536036) == "1140.5536" and num2str(27361.78648) == "27361.7865"
    assert num2str(0) == "0" and num2str(2.0) == "2" and num2str(0.5) == "0.5" and num2str(52) == "52"
    assert num2str(-1091961.0814, 6) == "-1.09196e+06"
    assert mat2str([1, 2.5]) == "[1 2.5]" and mat2str(3) == "3"
    assert [sp.num2ordinal(k) for k in (1, 2, 3, 4, 11, 12, 13, 21, 22, 101, 111)] == [
        "1st", "2nd", "3rd", "4th", "11th", "12th", "13th", "21st", "22nd", "101st", "111th"]  # fmt: skip
    with pytest.raises(SubEpiPredictError):
        sp.num2ordinal(0)


def test_panel_grid():
    expected = {1: (1, 1), 2: (2, 1), 3: (3, 1), 4: (2, 2), 6: (2, 3), 8: (4, 2), 9: (3, 3), 12: (4, 3), 16: (4, 4),
                24: (5, 5), 32: (6, 6)}  # fmt: skip
    for n, grid in expected.items():
        assert sp.subepidemicPanelGrid(n) == grid
    for n in range(1, 60):
        rows, cols = sp.subepidemicPanelGrid(n)
        assert rows * cols >= n
    assert factor(360) == [2, 2, 2, 3, 3, 5]
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicPanelGrid(0)


def test_parameter_label():
    assert sp.subepidemicParameterLabel("r_2", 24.23, 7.634, 102.4) == "r_2=24.2(95%CI:7.63,102)"
    assert sp.subepidemicParameterLabel("K_1", 45612.0, 39570.0, 52449.0) == "K_1=45600(95%CI:39600,52400)"
    assert sp.subepidemicParameterLabel("p_1", 0.46, 1.0, 0.0) == "p_1=0.46(95%CI:1,0)"
    assert sp.subepidemicParameterLabel("a", 1.2345e-5, 3.0e8, np.nan) == "a=1.23e-05(95%CI:3e+08,nan)"
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicParameterLabel("r", "x", 1, 2)


def test_estimate_line_keeps_the_vertical_limits():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.bar([1, 2, 3], [4, 9, 2])
    yl = ax.get_ylim()
    line = sp.subepidemicEstimateLine(2.5, ax)
    assert line is not None and ax.get_ylim() == yl and list(line.get_xdata()) == [2.5, 2.5]
    assert sp.subepidemicEstimateLine(np.nan, ax) is None


def test_struct_and_stream():
    s = sp.Struct(a=1)
    s.b = 2
    assert s["b"] == 2 and s.a == 1 and s.copy() == s
    with pytest.raises(AttributeError):
        s.missing
    sp.rng(4)
    a = sp.get_rng().random()
    sp.rng(4)
    assert sp.get_rng().random() == a
    own = np.random.default_rng(0)
    assert sp.get_rng(own) is own
