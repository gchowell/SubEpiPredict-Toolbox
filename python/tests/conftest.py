"""Shared set-up of the tests.

The tests import the package from the folder above this one, so they run from
a checkout without installing anything: ``cd python && python -m pytest``.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")  # no windows

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON_DIR = os.path.dirname(HERE)
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import numpy as np  # noqa: E402
import pytest  # noqa: E402

import subepipredict as sp  # noqa: E402


def pytest_configure(config):
    # The toolbox's own warnings (parameters at a bound, weak models ...) are
    # expected with the small settings of the tests; tests that check for one
    # use pytest.warns.
    config.addinivalue_line("filterwarnings", "ignore::subepipredict.errors.SubEpiPredictWarning")


@pytest.fixture(autouse=True)
def _clean_state():
    """Every test starts with seed 1, an empty warning log and no figures."""
    import matplotlib.pyplot as plt

    sp.rng(1)
    sp.subepidemicWarningLog("start")
    yield
    plt.close("all")


@pytest.fixture(scope="session")
def input_dir():
    folder = sp.example_input_dir()
    if folder is None:
        pytest.skip("The example datasets (ensemble n-subepidemic code v1.0/input) were not found.")
    return folder


@pytest.fixture(scope="session")
def covid_series(input_dir):
    """U.S. daily COVID-19 deaths to 05-11-2020 (column 52) as incidence: 75 values."""
    data = np.loadtxt(os.path.join(input_dir, "cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt"), ndmin=2)
    c = data[:, 51]
    return np.concatenate(([c[0]], np.diff(c)))


def identifier_of(excinfo):
    """Identifier of the SubEpiPredictError caught by pytest.raises."""
    return excinfo.value.identifier
