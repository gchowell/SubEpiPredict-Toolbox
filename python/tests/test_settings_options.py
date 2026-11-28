"""Fit settings, option objects, the checks made before a fit, the fitted series, the seed."""

import os

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError, SubEpiPredictWarning

# ---- fit settings -------------------------------------------------------------


def test_default_fit_settings():
    S = sp.subepidemicFitSettings()
    assert S.logScaleR and S.logScaleK and S.logScaleAlpha and S.countLikelihoodOnRawData
    assert S.rUpperBoundFactor == 10 and S.fitMaxAttempts == 3 and S.useParallel is False
    assert S.polish and S.polishGradientTol == 10 and S.polishBudget == 100 and S.agreementTol == 1e-2
    assert S.continuationStart and S.bootstrapStarts == 1 and S.rngSeed == 1
    assert S.maxSkippedRefitFraction == 0.05 and S.baselineWindow is None and S.PIIscale is None
    assert S.evidenceRatioLimit == 100 and S.firstThresholdStartFactor == 3 and S.notIdentifiedRatio == 0.5
    assert S.rLB is None and S.KUB is None and dict(S.optimOptions) == {}
    # the shipped options_fit_settings changes nothing
    assert sp.subepidemicFitSettings(sp.options_fit_settings) == S


def test_overrides_from_a_dict_or_a_function():
    S = sp.subepidemicFitSettings(dict(bootstrapStarts=3, rngSeed=None, KUB=5e6, useParallel=1))
    assert S.bootstrapStarts == 3 and S.rngSeed is None and S.KUB == 5e6 and S.useParallel is True
    S = sp.subepidemicFitSettings(lambda: {"polish": False})
    assert S.polish is False
    # a settings struct is passed through unchanged
    assert sp.subepidemicFitSettings(S).polish is False


@pytest.mark.parametrize(
    "override",
    [dict(polishTol=1), dict(polish=3), dict(agreementTol=0), dict(rUpperBoundFactor=-1), dict(fitMaxAttempts=0),
     dict(bootstrapStarts=1.5), dict(rLB=[1, np.nan]), dict(rngSeed=-3), dict(rngSeed=2.5), dict(baselineWindow=0),
     dict(maxSkippedRefitFraction=2), dict(firstThresholdStartFactor=0.5), dict(notIdentifiedRatio=0),
     dict(evidenceRatioLimit=1), dict(PIIscale=[1, 0, 1]), dict(PIIscale=[1, 0, 2, 0, 1, 1]), dict(optimOptions=3)],
)
def test_invalid_fit_settings_are_refused(override):
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicFitSettings(override)
    assert identifier_of(e) == "SubEpiPredict:FitSettings"


def test_fit_signature_ignores_report_only_settings():
    base = sp.subepidemicFitSignature()
    assert base.fitRevision == "m1" and "PIIscale" not in base.settings and "rngSeed" in base.settings
    same, differences = sp.subepidemicFitSignature(base, compare=True)
    assert same and differences == []
    same, _ = sp.subepidemicFitSignature(base, settings=dict(evidenceRatioLimit=1000, baselineWindow=14), compare=True)
    assert same
    same, differences = sp.subepidemicFitSignature(base, settings=dict(bootstrapStarts=3, KUB=5e6), compare=True)
    assert not same and differences == ["KUB", "bootstrapStarts"]
    same, differences = sp.subepidemicFitSignature(None, compare=True)
    assert not same and "no signature" in differences[0]
    old = sp.Struct(fitRevision="k", settings=base.settings)
    same, differences = sp.subepidemicFitSignature(old, compare=True)
    assert not same and differences[0].startswith("fit revision")


def test_fit_signature_survives_storage_as_text():
    from subepipredict.io import from_json, to_json

    signature = sp.subepidemicFitSignature(settings=dict(KUB=[2e6, 5e6], optimOptions={"FunctionTolerance": 1e-8}))
    back = from_json(to_json(signature))
    same, differences = sp.subepidemicFitSignature(
        back, settings=dict(KUB=[2e6, 5e6], optimOptions={"FunctionTolerance": 1e-8}), compare=True
    )
    assert same, differences


def test_seed_makes_the_stream_repeat():
    info = sp.subepidemicSeed("test", verbose=False)
    a = sp.get_rng().random(3)
    assert info.seed == 1 and info.state is not None
    sp.subepidemicSeed("test", verbose=False)
    assert np.array_equal(sp.get_rng().random(3), a)
    # no seed: the stream is left alone
    info = sp.subepidemicSeed("test", dict(rngSeed=None), verbose=False)
    assert info.seed is None and not np.array_equal(sp.get_rng().random(3), a)


# ---- options ------------------------------------------------------------------


def test_quick_start_options():
    o = sp.options()
    assert (o.cumulative1, o.outbreakx, o.caddate1, o.cadregion, o.caddisease, o.datatype) == (
        1, 52, "05-11-2020", "USA", "coronavirus", "deaths")  # fmt: skip
    assert (o.DT, o.datevecfirst1, o.datevecend1) == (1, [2020, 2, 27], [2022, 5, 9])
    assert (o.smoothfactor1, o.calibrationperiod1, o.method1, o.dist1) == (7, 90, 3, 3)
    assert (o.numstartpoints, o.B, o.npatches_fixed, o.topmodelsx, o.flag1, o.onset_fixed) == (20, 40, 2, 4, 1, 0)
    f = sp.options_forecast()
    assert (f.getperformance, f.deletetempfiles, f.forecastingperiod, f.weight_type1) == (1, 1, 4, -1)
    r = sp.options_Rt()
    assert (r.type_GId1, r.mean_GI1, r.var_GI1, r.minGIcoverage) == (1, 11.7, 3.7, 0.99)


def test_rules_of_options_m():
    o = sp.options()
    assert o.replace(method1=4).dist1 == 4 and o.replace(method1=1).dist1 == 1 and o.replace(method1=5).dist1 == 5
    assert o.replace(method1=0, dist1=2).dist1 == 2  # not set automatically for method 0
    assert o.replace(method1=6, dist1=6).dist1 == 6
    assert o.replace(npatches_fixed=1).topmodelsx == 1
    assert o.replace(onset_fixed=1).topmodelsx == 2 and o.replace(onset_fixed=1, npatches_fixed=6).topmodelsx == 4
    assert sp.options_wastewater().caddisease == "Wasterwater_viral_activity_level"
    assert sp.options_measles_SC().onset_fixed == 1 and sp.options_forecast_covid().forecastingperiod == 30


EXAMPLES = ["options", "options_covid", "options_tutorial", "options_SARSsingapore2003", "options_ebolaCongo1995",
            "options_ebolaUganda2022", "options_videoYoutube_Faded", "options_videoYoutube_Uptown_Funk", "options_mpx",
            "options_measles_SC", "options_acutelunginjury", "options_covid_wastewater", "options_wastewater",
            "options_copd", "options_diabetes", "options_yearly"]  # fmt: skip


# Two bundled example files begin with a column of dates written as text
# (9/21/25, 01/01/2022), so they are not numeric matrices.
DATE_COLUMN_FILES = {"options_measles_SC", "options_covid_wastewater"}


@pytest.mark.parametrize("name", EXAMPLES)
def test_every_example_passes_the_checks_and_finds_its_data(name, input_dir):
    o = getattr(sp, name)()
    sp.subepidemicCheckSettings(o)
    from subepipredict.utils import data_file_name, load_matrix

    filename = data_file_name(o.cumulative1, o.cadtemporal, o.caddisease, o.datatype, o.cadregion, o.caddate1)
    path = os.path.join(input_dir, filename)
    assert os.path.isfile(path), filename
    if name in DATE_COLUMN_FILES:
        with pytest.raises(SubEpiPredictError) as e:
            load_matrix(path)
        assert identifier_of(e) == "SubEpiPredict:DataFileNotNumeric" and filename in e.value.message
        return
    data = load_matrix(path)
    assert data.shape[1] >= o.outbreakx and data.shape[0] >= 10


def test_check_settings_lists_every_problem():
    bad = sp.options().as_dict()
    bad.update(method1=0, dist1=3, flag1=[1, 1], B=0, DT=2)
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicCheckSettings(bad)
    message = e.value.message
    assert identifier_of(e) == "SubEpiPredict:InvalidSettings" and "(4 problems)" in message
    assert "method1 = 0 with dist1 = 3 cannot run" in message and "flag1" in message and "B = 0" in message
    assert "DT = 2" in message


@pytest.mark.parametrize(
    "change, text",
    [(dict(method1=2), "Pearson"), (dict(method1=7), "method1 = 7 is not valid"), (dict(dist1=9), "dist1 = 9"),
     (dict(method1=1, dist1=3), "needs dist1 = 1"), (dict(method1=4, dist1=3), "needs dist1 = 4"),
     (dict(method1=6, dist1=0), "Laplace"), (dict(npatches_fixed=0), "npatches_fixed"), (dict(onset_fixed=2), "onset_fixed"),
     (dict(smoothfactor1=0), "smoothfactor1"), (dict(calibrationperiod1=1), "calibrationperiod1"),
     (dict(numstartpoints=0), "numstartpoints"), (dict(topmodelsx=0), "topmodelsx"), (dict(cumulative1=2), "cumulative1"),
     (dict(outbreakx=[1, 0]), "outbreakx")],
)
def test_check_settings_single_problems(change, text):
    S = sp.options().as_dict()
    S.update(change)
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicCheckSettings(S)
    assert text in e.value.message and e.value.message.startswith("The settings cannot run. ")


def test_check_settings_accepts_valid_combinations_and_reports_missing_fields():
    S = sp.options().as_dict()
    for method1, dist1 in ((0, 0), (0, 1), (0, 2), (1, 1), (3, 3), (4, 4), (5, 5), (6, 6), (6, 1)):
        S.update(method1=method1, dist1=dist1)
        sp.subepidemicCheckSettings(S)
    del S["B"]
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicCheckSettings(S)
    assert "B is missing." in e.value.message


# ---- fitted series ------------------------------------------------------------

SERIES = np.array([1.0, 0, 3, 5, 9, 14, 12, 20, 18, 11])


def test_count_likelihoods_are_fitted_to_the_raw_counts():
    for method1 in (1, 3, 4, 5):
        y, info = sp.subepidemicFitSeries(SERIES, 7, method1)
        assert np.array_equal(y, SERIES) and info.usedRawData and info.smoothingSpan == 1
    y, info = sp.subepidemicFitSeries(SERIES, 7, 3, dict(countLikelihoodOnRawData=False))
    assert info.smoothedCountLikelihood and np.allclose(y, sp.movmean(SERIES, 7))


def test_least_squares_is_fitted_to_the_moving_average():
    for method1 in (0, 6):
        y, info = sp.subepidemicFitSeries(SERIES, 3, method1)
        assert np.allclose(y, sp.movmean(SERIES, 3)) and info.smoothingSpan == 3 and not info.usedRawData
    y, info = sp.subepidemicFitSeries(SERIES, 1, 0)
    assert np.array_equal(y, SERIES) and info.usedRawData


def test_fit_series_errors_and_warning():
    negative = SERIES.copy()
    negative[3] = -2
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicFitSeries(negative, 7, 3)
    assert identifier_of(e) == "SubEpiPredict:NegativeCounts"
    assert sp.subepidemicFitSeries(negative, 1, 0)[0][3] == -2  # least squares accepts it
    with pytest.warns(SubEpiPredictWarning, match="NonIntegerCounts"):
        _, info = sp.subepidemicFitSeries(SERIES + 0.5, 7, 1)
    assert info.nNonInteger == SERIES.size
    for args, identifier in (
        (([], 7, 3), "SubEpiPredict:InvalidData"),
        ((SERIES, 7, 9), "SubEpiPredict:InvalidMethod"),
        ((SERIES, 0.5, 3), "SubEpiPredict:InvalidSmoothing"),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicFitSeries(*args)
        assert identifier_of(e) == identifier


# ---- warnings -----------------------------------------------------------------


def test_warning_log_counts_kinds_and_prints_a_summary(capsys):
    sp.subepidemicWarningLog("start")
    with pytest.warns(SubEpiPredictWarning) as record:
        sp.subepidemicWarning("SubEpiPredict:TestA", "first %d", 1)
        sp.subepidemicWarning("SubEpiPredict:TestA", "second %d", 2)
        sp.subepidemicWarning("SubEpiPredict:TestB", "other")
    assert record[0].message.identifier == "SubEpiPredict:TestA"
    log = sp.subepidemicWarningLog("get")
    assert [(e["id"], e["count"], e["message"]) for e in log] == [
        ("SubEpiPredict:TestA", 2, "first 1"), ("SubEpiPredict:TestB", 1, "other")]  # fmt: skip
    sp.subepidemicWarningLog("summary", "my run")
    out = capsys.readouterr().out
    assert "my run finished with 2 kind(s) of warning" in out and "(2 times; first:) first 1" in out
    sp.subepidemicWarningLog("start")
    assert sp.subepidemicWarningLog("get") == []
    with pytest.raises(SubEpiPredictError):
        sp.subepidemicWarningLog("nothing")
