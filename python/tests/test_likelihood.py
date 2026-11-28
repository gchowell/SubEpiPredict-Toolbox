"""Negative-binomial likelihood and the fitting objective."""

import math

import numpy as np
import pytest
from conftest import identifier_of
from scipy.special import gammaln
from scipy.stats import nbinom, poisson

import subepipredict as sp
from subepipredict import SubEpiPredictError

RNG = np.random.default_rng(7)
Y = RNG.integers(0, 3000, 60).astype(float)
MU = RNG.uniform(0.01, 3000, 60)


def reference_nll(y, mu, alpha, d, method):
    var = {3: mu + alpha * mu, 4: mu + alpha * mu**2, 5: mu + alpha * mu**d}[method]
    p = mu / var
    r = mu * p / (1 - p)
    return -(nbinom.logpmf(y, r, p).sum() + gammaln(y + 1).sum())


@pytest.mark.parametrize(
    "method, alpha, d",
    [(3, 0.5, 1), (3, 30.0, 1), (3, 1e5, 1), (4, 0.02, 1), (4, 50, 1), (5, 0.3, 1.4), (5, 2, 3), (5, 0.4, 0.5)],
)
def test_matches_the_negative_binomial_of_scipy(method, alpha, d):
    value = sp.subepidemicNegativeBinomialNLL(Y, MU, alpha, d, method)
    assert value == pytest.approx(reference_nll(Y, MU, alpha, d, method), rel=1e-11)


def test_small_dispersion_approaches_the_poisson_limit():
    poissonValue = -(poisson.logpmf(Y, MU).sum() + gammaln(Y + 1).sum())
    assert sp.subepidemicNegativeBinomialNLL(Y, MU, 0, 1, 3) == pytest.approx(poissonValue, rel=1e-13)
    assert sp.subepidemicNegativeBinomialNLL(Y, MU, 1e-10, 1, 4) == pytest.approx(poissonValue, rel=1e-6)


def test_row_and_column_inputs_and_large_counts():
    a = sp.subepidemicNegativeBinomialNLL(Y.reshape(-1, 1), MU.reshape(1, -1), 0.7, 1, 3)
    assert a == sp.subepidemicNegativeBinomialNLL(Y, MU, 0.7, 1, 3)
    big = sp.subepidemicNegativeBinomialNLL([1e7, 3e8], [1.1e7, 2.9e8], 1e-3, 1, 3)
    assert math.isfinite(big)


def test_fractional_data_keep_the_floor_product_convention():
    y, mu, alpha = np.array([2.6]), np.array([3.0]), 0.5
    r = mu / alpha
    beta = alpha
    expected = -(gammaln(r + 2) - gammaln(r) + y * np.log(beta / (1 + beta)) - r * np.log(1 + beta))[0]
    assert sp.subepidemicNegativeBinomialNLL(y, mu, alpha, 1, 3) == pytest.approx(expected, rel=1e-12)
    # at alpha = 0 fractional data have no finite legacy objective
    assert sp.subepidemicNegativeBinomialNLL(y, mu, 0, 1, 3) == math.inf


@pytest.mark.parametrize(
    "y, mu, alpha, d",
    [([1, 2], [0, 2], 0.5, 1), ([1, -2], [1, 2], 0.5, 1), ([1, 2], [1, np.nan], 0.5, 1), ([1, 2], [1, 2], -0.1, 1),
     ([1, 2], [1, 2], np.inf, 1), ([1, 2], [1, 2], 0.5, np.nan)],
)
def test_invalid_numbers_give_infinity(y, mu, alpha, d):
    assert sp.subepidemicNegativeBinomialNLL(y, mu, alpha, d, 5) == math.inf


def test_shape_and_method_errors():
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicNegativeBinomialNLL([1, 2, 3], [1, 2], 0.5, 1, 3)
    assert identifier_of(e) == "SubEpiPredict:NBInputShape"
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicNegativeBinomialNLL([1, 2], [1, 2], 0.5, 1, 2)
    assert identifier_of(e) == "SubEpiPredict:NBMethod"


# ---- objective ----------------------------------------------------------------

T = np.arange(0, 30, dtype=float)
Z = np.array([0.4, 0.9, 1.0, 2000.0, 0.8, 1.2])  # one sub-epidemic: r p a K alpha d
CURVE = sp.simulateSubepidemic(T, 3, 1, 0, 0, 1, [0.4], [0.9], [1.0], [2000.0])[2]
DATA = np.round(CURVE * (1 + 0.2 * np.sin(T)))


def context(method1, ydata=DATA, **changes):
    c = sp.Struct(timevect=T, ydata=ydata, I0=3.0, npatches=1, onset_fixed=0, onset_thr=0.0, flag1=1, method1=method1,
                  odeOptions=sp.subepidemicOdeOptions())  # fmt: skip
    c.update(changes)
    return c


def test_objective_values_of_every_method():
    v, detail = sp.subepidemicObjective(Z, context(0))
    assert detail.success and v == pytest.approx(np.sum((DATA - CURVE) ** 2), rel=1e-10)
    assert np.allclose(detail.curve, CURVE) and detail.nActive == 1
    v, _ = sp.subepidemicObjective(Z, context(1))
    assert v == pytest.approx(-np.sum(DATA * np.log(CURVE) - CURVE), rel=1e-10)
    v, _ = sp.subepidemicObjective(Z, context(6))
    assert v == pytest.approx(np.sum(np.abs(DATA - CURVE)), rel=1e-10)
    for method1 in (3, 4, 5):
        v, _ = sp.subepidemicObjective(Z, context(method1))
        assert v == pytest.approx(sp.subepidemicNegativeBinomialNLL(DATA, CURVE, Z[-2], Z[-1], method1), rel=1e-12)


def test_objective_marks_invalid_vectors_with_the_penalty():
    v, detail = sp.subepidemicObjective(Z[:-1], context(0))
    assert v == 1e10 and not detail.success and "wrong size" in detail.message
    bad = Z.copy()
    bad[0] = np.nan
    assert sp.subepidemicObjective(bad, context(0))[0] == 1e10
    # a trajectory that cannot be integrated (finite-time blow-up of the GGM with p > 1)
    v, detail = sp.subepidemicObjective([5.0, 3.0, 1, 20, 0, 0], context(0, flag1=0, timevect=np.arange(0, 400.0),
                                                                        ydata=np.ones(400)))  # fmt: skip
    assert v == 1e10 and not detail.success and detail.message.startswith("SubEpiPredict:")
    # unknown estimation method
    v, detail = sp.subepidemicObjective(Z, context(2))
    assert v == 1e10 and "InvalidMethod" in detail.message


def test_objective_rejects_an_identically_zero_trajectory():
    v, detail = sp.subepidemicObjective([0.0, 1, 1, 100, 0, 0], context(0, I0=0.0))
    assert v == 1e10 and not detail.success and "identically zero" in detail.message


def test_a_valid_poor_fit_is_not_penalized():
    v, detail = sp.subepidemicObjective(Z, context(0, ydata=DATA * 1e6))
    assert detail.success and v > 1e10
