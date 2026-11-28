"""AICc, the number of parameters, and the error structures."""

import math

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError, SubEpiPredictWarning

# ---- number of parameters and AICc --------------------------------------------


@pytest.mark.parametrize("flag1, per", [(0, 2), (1, 3), (2, 4), (3, 2), (4, 3), (5, 2)])
def test_parameters_per_sub_epidemic(flag1, per):
    assert sp.get_nparams(1, 1, 1, flag1, 1, 0) == per
    assert sp.get_nparams(1, 1, 3, flag1, 1, 1) == 3 * per  # synchronous: no threshold
    assert sp.get_nparams(1, 1, 3, flag1, 1, 0) == 3 * per + 1  # asynchronous: one threshold
    assert sp.get_nparams(1, 1, 1, flag1, 0, 0) == per + 1  # initial value estimated


@pytest.mark.parametrize("method1, dist1, extra", [(0, 0, 1), (0, 1, 0), (0, 2, 0), (1, 1, 0), (3, 3, 1), (4, 4, 1),
                                                   (5, 5, 2), (6, 6, 1)])  # fmt: skip
def test_parameters_of_the_error_structure(method1, dist1, extra):
    assert sp.get_nparams(method1, dist1, 2, 1, 1, 0) == 7 + extra


def test_aicc_formulas():
    n, k = 75, 8  # two GLM sub-epidemics, asynchronous, method 3
    AICc, part1, part2, numparams = sp.getAICc(3, 3, 2, 1, 1, -545989.6316, n, 0)
    assert numparams == k
    assert part1 == pytest.approx(2 * -545989.6316) and part2 == pytest.approx(2 * k + 2 * k * (k + 1) / (n - k - 1))
    assert AICc == pytest.approx(part1 + part2)
    # least squares: n*log(SSE)
    AICc, part1, part2, numparams = sp.getAICc(0, 0, 1, 3, 1, 1234.5, 30, 0)
    assert numparams == 3 and part1 == pytest.approx(30 * math.log(1234.5))
    assert AICc == pytest.approx(30 * math.log(1234.5) + 6 + 24 / 26)
    # least absolute deviations: 2n*log(2*SAD/n) + 2n
    AICc, part1, part2, _ = sp.getAICc(6, 6, 1, 3, 1, 300.0, 30, 0)
    assert part1 == pytest.approx(60 * math.log(20.0)) and part2 == pytest.approx(60 + 6 + 24 / 26)


def test_aicc_is_infinite_when_the_data_cannot_support_the_model():
    # 3 GLM sub-epidemics with least squares: 11 parameters
    assert sp.get_nparams(0, 0, 3, 1, 1, 0) == 11
    for n in (10, 12):
        AICc, part1, part2, _ = sp.getAICc(0, 0, 3, 1, 1, 50.0, n, 0)
        assert AICc == math.inf and math.isnan(part1) and part2 == math.inf
    assert math.isfinite(sp.getAICc(0, 0, 3, 1, 1, 50.0, 13, 0)[0])


# ---- error structures -----------------------------------------------------------

CUM = np.cumsum(np.array([5.0, 20, 60, 150, 300, 420, 380, 200]))
MEAN = np.diff(CUM)
M = 40000


def draws(dist1, factor1=None, d=None, seed=11):
    return sp.AddPoissonError(CUM, M, dist1, factor1, d, rng=np.random.default_rng(seed))


def test_first_row_is_not_perturbed_and_shape():
    c = draws(1)
    assert c.shape == (CUM.size, M) and np.all(c[0] == CUM[0])
    assert sp.AddPoissonError(CUM, 0, 1).shape == (CUM.size, 0)
    assert sp.AddPoissonError([7.0], 3, 3, 0.5).tolist() == [[7.0, 7.0, 7.0]]


@pytest.mark.parametrize(
    "dist1, factor1, d, variance",
    [
        (0, 12.0, None, lambda m: np.full(m.shape, 144.0)),
        (1, None, None, lambda m: m),
        (2, 4.0, None, lambda m: 4.0 * m),
        (3, 2.5, None, lambda m: m + 2.5 * m),
        (4, 0.05, None, lambda m: m + 0.05 * m**2),
        (5, 0.3, 1.5, lambda m: m + 0.3 * m**1.5),
        (6, 9.0, None, lambda m: np.full(m.shape, 2 * 81.0)),
    ],
)
def test_mean_and_variance_of_every_error_structure(dist1, factor1, d, variance):
    c = draws(dist1, factor1, d)[1:]
    expected = variance(MEAN)
    se = np.sqrt(expected / M)
    assert np.all(np.abs(c.mean(axis=1) - MEAN) < 6 * se)
    assert np.allclose(c.var(axis=1, ddof=1), expected, rtol=0.08)
    if dist1 in (1, 2, 3, 4, 5):
        assert np.all(c >= 0) and np.all(c == np.round(c))


def test_zero_mean_and_zero_dispersion():
    flat = np.array([3.0, 3.0, 10.0, 10.0])
    for dist1, factor1 in ((1, None), (2, 3.0), (3, 2.0), (4, 0.1), (5, 0.1)):
        c = sp.AddPoissonError(flat, 50, dist1, factor1, 1.3)
        assert np.all(c[1] == 0) and np.all(c[3] == 0) and np.all(np.isfinite(c))
    # zero overdispersion is the Poisson limit
    for dist1, factor1 in ((2, 1.0), (3, 0.0), (4, 0.0), (5, 0.0)):
        c = sp.AddPoissonError(CUM, 20000, dist1, factor1, 1.3, rng=np.random.default_rng(2))[1:]
        assert np.allclose(c.var(axis=1, ddof=1), MEAN, rtol=0.1)


def test_same_seed_same_curves():
    sp.rng(5)
    a = sp.AddPoissonError(CUM, 7, 3, 1.5)
    sp.rng(5)
    assert np.array_equal(sp.AddPoissonError(CUM, 7, 3, 1.5), a)
    assert not np.array_equal(sp.AddPoissonError(CUM, 7, 3, 1.5), a)


def test_warnings_for_underdispersion_and_decreasing_curves():
    with pytest.warns(SubEpiPredictWarning, match="UnderdispersedFactor"):
        c = sp.AddPoissonError(CUM, 5000, 2, 0.4, rng=np.random.default_rng(1))
    assert np.allclose(c[1:].var(axis=1, ddof=1), MEAN, rtol=0.2)  # Poisson instead
    down = np.array([10.0, 30.0, 25.0, 40.0])
    with pytest.warns(SubEpiPredictWarning, match="DecreasingCurve"):
        c = sp.AddPoissonError(down, 100, 1)
    assert np.all(c[2] == 0)
    # ODE noise on a plateau is passed over silently
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sp.AddPoissonError(np.array([10.0, 1000.0, 1000.0 - 1e-7]), 10, 3, 1.0)


@pytest.mark.parametrize(
    "args",
    [([], 1, 1), ([1, np.nan], 1, 1), ([1, 2], -1, 1), ([1, 2], 1.5, 1), ([1, 2], 1, 7), ([1, 2], 1, 3),
     ([1, 2], 1, 3, -0.5), ([1, 2], 1, 0, np.nan), ([1, 2], 1, 5, 0.5, 0), ([1, 2], 1, 3, [1, 2])],
)
def test_error_generator_checks_its_inputs(args):
    with pytest.raises(SubEpiPredictError) as e:
        sp.AddPoissonError(*args)
    assert identifier_of(e) == "SubEpiPredict:ErrorGenerator"


def test_mean_variance_ratio():
    data = np.array([2.0, 4, 6, 8, 10, 12, 14, 16, 18])
    ratios, bins = sp.getMeanVarianceRatio(data, 3, 2)  # moving windows starting at every point
    assert ratios.shape == (6, 3) and bins.tolist() == [1, 2, 3, 4, 5, 6]
    assert ratios[0].tolist() == [1.0, 4.0, 4.0]  # var([2 4 6]) = 4, mean 4
    ratios, bins = sp.getMeanVarianceRatio(data, 3, 1)  # consecutive windows
    assert ratios[:, 1].tolist() == [4.0, 10.0] and bins.tolist() == [2.0, 5.0]
