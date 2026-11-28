"""End-to-end: fit, plots, forecast, ensembles, reproduction number, rolling forecasts.

Small settings keep the runs short: the bundled U.S. COVID-19 deaths series,
logistic sub-epidemics, least squares, few start points and bootstrap refits.
"""

import glob
import os

import numpy as np
import pytest
from conftest import identifier_of

import subepipredict as sp
from subepipredict import SubEpiPredictError, SubEpiPredictWarning
from subepipredict.io import OutputNames, readtable

AREA = 52
DATE = "04-20-2020"  # 54 data rows; the last 30 are used


# The default lower bound of r is the mean of the first two observations divided
# by 200. A window that starts in the middle of an epidemic then starts above
# the growth rate of a logistic curve (p fixed at 1), so the tests set a floor.
FS = dict(rLB=0.01)


def small_options(**changes):
    base = dict(method1=0, dist1=0, flag1=sp.LM, smoothfactor1=1, calibrationperiod1=30, numstartpoints=2, B=6,
                npatches_fixed=2, topmodelsx=3, caddate1=DATE)  # fmt: skip
    base.update(changes)
    return sp.options().replace(**base)


def one_model_options(**changes):
    return small_options(**{**dict(npatches_fixed=1, numstartpoints=3, B=8), **changes})


def files(folder, pattern):
    return sorted(glob.glob(os.path.join(glob.escape(folder), pattern)))


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    """One fit with up to two sub-epidemics, shared by the tests below."""
    out = str(tmp_path_factory.mktemp("fit"))
    opts = small_options()
    results = sp.Run_Fit_subepidemicFramework(AREA, DATE, options=opts, fit_settings=FS, output_dir=out, plot=False,
                                              verbose=False)  # fmt: skip
    return opts, out, results


# ---- fit ----------------------------------------------------------------------------


@pytest.mark.slow
def test_fit_returns_and_saves_the_ranked_models(fitted):
    opts, out, results = fitted
    assert len(results) == 1
    r = results[0]
    assert r.outbreakx == AREA and r.topmodelsx == 3 and len(r.ranked) == 3
    assert np.all(np.diff(r.RMSES[:, 2]) >= 0) and r.relativelik_i[0] == 1
    assert r.evidenceRatio[0] == 1 and np.all(np.diff(r.evidenceRatio) >= 0)
    assert r.data1.shape == (30, 2) and r.data1[0, 0] == 0 and r.timelags == 24
    names = OutputNames(opts, AREA, DATE, out)
    assert os.path.basename(names.abc_mat()) == (
        "ABC-ensem-npatchesfixed-2-onsetfixed-0-smoothing-1-daily-coronavirus-deaths-USA-state-52-04-20-2020"
        "-flag1-3-method-0-dist-0-calibrationperiod-30.mat")  # fmt: skip
    assert os.path.isfile(names.abc_mat())
    for k, ranked in enumerate(r.ranked, start=1):
        assert ranked.file == names.ranked_mat(k) and ranked.file.endswith("-calibrationperiod-30-rank-%d.mat" % k)
        saved = sp.load_mat(ranked.file)
        n = int(saved.npatches[0, 0])
        assert saved.Phatss.shape == (6, 4 * n + 2) and saved.curves.shape == (30, 6) and saved.data1.shape == (30, 2)
        assert np.allclose(saved.P0.ravel(), ranked.P0) and float(saved.AICc_best[0, 0]) == r.RMSES[k - 1, 2]
        assert int(saved.timelags[0, 0]) == 24 and int(saved.M[0, 0]) == 6 and saved.caddate1 == DATE
        assert np.allclose(saved.bootstrapDiagnostics.searchLB.ravel(), ranked.bootstrapDiagnostics.searchLB)
        assert saved.fitSignature.fitRevision == "m1"
    assert not files(out, "*.csv")  # the fit writes MAT files only


@pytest.mark.slow
def test_fit_is_repeatable_and_several_areas_are_fitted_in_turn(tmp_path):
    opts = one_model_options()
    kw = dict(options=opts, fit_settings=FS, plot=False, verbose=False)
    a = sp.Run_Fit_subepidemicFramework(AREA, DATE, output_dir=str(tmp_path / "a"), **kw)
    b = sp.Run_Fit_subepidemicFramework([33, AREA], DATE, output_dir=str(tmp_path / "b"), **kw)
    assert [r.outbreakx for r in b] == [33, AREA]
    # the seed is set per area: an area fitted after another equals the same area fitted alone
    assert np.array_equal(a[0].RMSES, b[1].RMSES) and np.array_equal(a[0].ranked[0].Phatss, b[1].ranked[0].Phatss)
    assert np.array_equal(a[0].ranked[0].curves, b[1].ranked[0].curves)
    assert len(files(str(tmp_path / "b"), "modifiedLogisticPatch-*-state-33-*-rank-1.mat")) == 1
    # without a seed, runs differ
    c = sp.Run_Fit_subepidemicFramework(AREA, DATE, options=opts, fit_settings=dict(rngSeed=None, **FS),
                                        output_dir=str(tmp_path / "c"), plot=False, verbose=False)  # fmt: skip
    assert not np.array_equal(a[0].ranked[0].curves, c[0].ranked[0].curves)


def test_fit_stops_before_fitting_when_the_settings_cannot_run(tmp_path):
    out = str(tmp_path)
    with pytest.raises(SubEpiPredictError) as e:
        sp.Run_Fit_subepidemicFramework(AREA, DATE, options=small_options(method1=0, dist1=3), output_dir=out)
    assert identifier_of(e) == "SubEpiPredict:InvalidSettings"
    with pytest.raises(SubEpiPredictError) as e:
        sp.Run_Fit_subepidemicFramework(AREA, DATE, options=small_options(calibrationperiod1=5), output_dir=out)
    assert identifier_of(e) == "SubEpiPredict:TooManyParameters"
    with pytest.raises(SubEpiPredictError) as e:
        sp.Run_Fit_subepidemicFramework(99, DATE, options=small_options(), output_dir=out, verbose=False)
    assert identifier_of(e) == "SubEpiPredict:InvalidSettings" and "52 column(s)" in e.value.message
    with pytest.raises(FileNotFoundError):
        sp.Run_Fit_subepidemicFramework(AREA, "01-01-1999", options=small_options(), output_dir=out, verbose=False)
    assert not os.listdir(out)


# ---- plots of the fit ------------------------------------------------------------------


@pytest.mark.slow
def test_plot_fit_writes_the_parameter_tables(fitted):
    import matplotlib.pyplot as plt

    opts, out, results = fitted
    performance = sp.plotFit_subepidemicFramework(AREA, DATE, options=opts, output_dir=out, verbose=False)
    assert performance.shape == (3, 7) and performance[:, 0].tolist() == [1, 2, 3]
    assert np.allclose(performance[:, 5], results[0].RMSES[:3, 2]) and performance[0, 6] == 1
    assert np.all((performance[:, 3] >= 0) & (performance[:, 3] <= 100)) and np.all(performance[:, 4] > 0)
    assert set(plt.get_fignums()) == {301, 302, 303, 304, 400}
    assert len(plt.figure(304).axes) == 9  # three ranked models: fit, components, residuals
    tail = "-daily-coronavirus-deaths-USA-area-52-04-20-2020.csv"
    for k in (1, 2, 3):
        n = results[0].ranked[k - 1].npatches
        for letter in "rpaK":
            names_, table = readtable(os.path.join(
                out, "param-%s-ranked(%d)-onsetfixed-0-flag1-3-method-0-dist-0-horizon-0%s" % (letter, k, tail)))  # fmt: skip
            assert names_ == ["sub-epidemic", letter + " mean", letter + " 95% CI LB", letter + " 95% CI UB",
                              letter + " SCI", letter + " estimate", letter + " median"]  # fmt: skip
            assert table.shape == (n, 7)
        P0 = results[0].ranked[k - 1].P0
        assert np.allclose(table[:, 5], P0[3 * n : 4 * n])  # K estimate = the ranked fit
        assert np.all(table[:, 2] <= table[:, 6]) and np.all(table[:, 6] <= table[:, 3])
        assert os.path.isfile(os.path.join(out, "PII-ranked(%d)-onsetfixed-0-flag1-3-method-0-dist-0-horizon-0%s" % (k, tail)))
        assert os.path.isfile(os.path.join(out, "doublingTimes-ranked(%d)-onsetfixed-0-flag1-3-method-0-dist-0-horizon-0%s" % (k, tail)))
    assert not files(out, "param-NB-*")  # least squares has no alpha
    names_, table = readtable(os.path.join(out, "performance-calibration-topRanked-onsetfixed-0-flag1-3-method-0-dist-0" + tail))
    assert names_[0] == "i_th-ranked model" and np.allclose(table, performance)
    lines = open(os.path.join(out, "PII-topRanked-onsetfixed-0-flag1-3-method-0-dist-0" + tail)).read().splitlines()
    assert lines[0].startswith("ranked_model,n_subepidemics,parameter") and lines[1].startswith("1,")


@pytest.mark.slow
def test_plot_rankings(fitted):
    import matplotlib.pyplot as plt

    opts, out, results = fitted
    top = sp.plotRankings_subepidemicFramework(AREA, DATE, options=opts, output_dir=out, verbose=False)
    assert np.allclose(top[:, 1], results[0].RMSES[:3, 2]) and np.allclose(top[:, 2], results[0].relativelik_i[:3])
    assert set(plt.get_fignums()) == {99, 100} and len(plt.figure(100).axes) == 3
    names_, table = readtable(os.path.join(
        out, "AICc-topRanked-onsetfixed-0-flag1-3-method-0-dist-0-daily-coronavirus-deaths-USA-area-52-04-20-2020.csv"))  # fmt: skip
    assert names_ == ["i_th-ranked model", "AICc", "RelativeLikelihood"] and np.allclose(table, top)
    assert len(files(out, "AICc_bests-*-dateini-04-20-20-*-topmodels-3.mat")) == 1
    # asking for more top models than the screen kept, or than were bootstrapped
    many = opts.replace(topmodelsx=500)
    with pytest.warns(SubEpiPredictWarning, match="TopModelsCapped"):
        top = sp.plotRankings_subepidemicFramework(AREA, DATE, options=many, output_dir=out, plot=False, verbose=False)
    assert top.shape[0] == results[0].RMSES.shape[0]
    with pytest.warns(SubEpiPredictWarning, match="TopModelsCapped"):
        performance = sp.plotFit_subepidemicFramework(AREA, DATE, options=many, output_dir=out, plot=False, verbose=False)
    assert performance.shape[0] == 3


# ---- forecasts ----------------------------------------------------------------------------


@pytest.mark.slow
def test_forecast_scores_ranked_models_ensembles_and_baseline(fitted):
    import matplotlib.pyplot as plt

    opts, out, results = fitted
    ranked, ensembles, baseline, byHorizon = sp.plotForecast_subepidemicFramework(
        AREA, DATE, 5, -1, options=opts, output_dir=out, verbose=False
    )
    assert ranked.shape == (3, 7) and ranked[:, 0].tolist() == [1, 2, 3]
    assert np.allclose(ranked[:, 5], results[0].RMSES[:3, 2])
    assert ensembles.shape == (2, 5) and ensembles[:, 0].tolist() == [2, 3]
    assert baseline.shape == (4,) and np.all(np.isfinite(baseline))
    # one row per model and horizon: baseline 0, ranked 1-3, ensembles 102 and 103
    assert byHorizon.shape == (6 * 5, 7)
    assert sorted(set(byHorizon[:, 0].tolist())) == [0, 1, 2, 3, 102, 103]
    assert np.all(byHorizon[byHorizon[:, 0] == 0, 6] == 1)  # the baseline relative to itself
    last = byHorizon[byHorizon[:, 1] == 5]
    assert np.allclose(last[last[:, 0] == 0, 2:6].ravel(), baseline)
    assert np.allclose(last[(last[:, 0] >= 1) & (last[:, 0] <= 3), 2], ranked[:, 1])
    assert np.allclose(last[last[:, 0] > 100, 5], ensembles[:, 4])
    assert np.allclose(byHorizon[:, 6], byHorizon[:, 5] / np.tile(byHorizon[byHorizon[:, 0] == 0, 5], 6))
    # the baseline is the mean of the last seven observations
    data1 = results[0].data1
    later = np.diff(np.loadtxt(os.path.join(sp.example_input_dir(), "cumulative-daily-coronavirus-deaths-USA-05-09-2022.txt"))[:, 51])[53:58]
    assert baseline[0] == pytest.approx(np.mean(np.abs(later - data1[-7:, 1].mean())))
    assert {10, 11, 500, 400, 401} <= set(plt.get_fignums()) and len(plt.figure(500).axes) == 2

    tail = "-daily-coronavirus-deaths-USA-area-52-04-20-2020.csv"
    head = "-onsetfixed-0-flag1-3-method-0-dist-0-horizon-5"
    names_, table = readtable(os.path.join(out, "ranked(1)" + head + tail))
    assert names_ == ["year", "month", "day", "data", "median", "LB", "UB"] and table.shape == (35, 7)
    assert table[0, :3].tolist() == [2020, 3, 22] and table[-1, :3].tolist() == [2020, 4, 25]
    assert np.allclose(table[:30, 3], data1[:, 1]) and np.allclose(table[30:, 3], later)
    assert np.all(table[:, 5] <= table[:, 6]) and np.all(table[:, 5] >= 0)  # limits below zero are set to 0
    names_, table = readtable(os.path.join(out, "Ensemble(3)" + head + "-weighttype--1" + tail))
    assert table.shape == (35, 7) and np.allclose(table[30:, 3], later)
    names_, q = readtable(os.path.join(out, "quantile-ranked(2)" + head + tail))
    assert names_[0] == "Q_0.010" and names_[-1] == "Q_0.990" and q.shape == (35, 23) and np.all(np.diff(q, axis=1) >= 0)
    assert readtable(os.path.join(out, "quantileTimes-Ensemble(2)" + head + "-weighttype--1" + tail))[1].shape == (35, 23)
    for name in ("performance-forecasting-topRanked" + head, "performance-forecasting-Ensemble" + head + "-weight_type--1",
                 "performance-forecasting-byHorizon" + head + "-weight_type--1",
                 "performance-calibration-Ensemble" + head + "-weight_type--1",
                 "performance-calibration-topRanked-onsetfixed-0-flag1-3-method-0-dist-0",
                 "doublingTimes-ranked(1)" + head, "doublingTimes-Ensemble(3)" + head + "-weighttype--1"):  # fmt: skip
        assert os.path.isfile(os.path.join(out, name + tail)), name
    assert not files(out, "Forecast-*.mat")  # temporary files are deleted


@pytest.mark.slow
def test_forecast_is_repeatable_and_weights_change_only_the_ensembles(fitted):
    opts, out, _ = fitted
    kw = dict(options=opts, output_dir=out, plot=False, verbose=False)
    r1, e1, b1, h1 = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, -1, **kw)
    r2, e2, b2, h2 = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, -1, **kw)
    assert np.array_equal(h1, h2)  # the seed is reset at the start of the function
    ra, ea, _, _ = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 1, **kw)
    assert np.array_equal(ra, r1) and not np.array_equal(ea, e1)
    with pytest.warns(SubEpiPredictWarning, match="WeightCodeZero"):
        _, e0, _, _ = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 0, **kw)
    assert np.array_equal(e0, ea)  # code 0 is code 1
    _, ew, _, _ = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 2, **kw)
    _, e3, _, _ = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 3, **kw)
    assert np.array_equal(ew, e3) and not np.array_equal(ew, e1)
    assert len(files(out, "Ensemble(2)-*-horizon-5-weighttype-2-*.csv")) == 1
    with pytest.raises(SubEpiPredictError) as e:
        sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 7, **kw)
    assert identifier_of(e) == "SubEpiPredict:EnsembleWeights"
    # a horizon beyond the data that exist
    with pytest.raises(SubEpiPredictError) as e:
        sp.plotForecast_subepidemicFramework(AREA, DATE, 2000, -1, **kw)
    assert identifier_of(e) == "SubEpiPredict:ForecastTooLong"


@pytest.mark.slow
def test_forecast_without_scores_and_with_kept_curves(fitted):
    opts, out, _ = fitted
    fopts = sp.options_forecast().replace(getperformance=0, deletetempfiles=0, forecastingperiod=7)
    ranked, ensembles, baseline, byHorizon = sp.plotForecast_subepidemicFramework(
        AREA, DATE, options=opts, options_forecast=fopts, output_dir=out, plot=False, verbose=False
    )
    assert ranked.shape == (0, 7) and ensembles.shape == (0, 5) and baseline.size == 0 and byHorizon.shape == (0, 7)
    # without scores the name of the ranked forecast file has no horizon
    path = os.path.join(out, "ranked(1)-onsetfixed-0-flag1-3-method-0-dist-0-daily-coronavirus-deaths-USA-area-52-04-20-2020.csv")
    names_, table = readtable(path)
    assert table.shape == (37, 7) and np.all(np.isnan(table[30:, 3])) and np.all(np.isfinite(table[:, 4:]))
    kept = files(out, "Forecast-modifiedLogisticPatch-*-forecastingperiod-7-rank-*.mat")
    assert len(kept) == 3
    saved = sp.load_mat(kept[0])
    assert saved.curvesforecasts1.shape == (37, 6) and saved.curvesforecasts2.shape == (37, 120)
    assert np.all(saved.curvesforecasts2[0] == saved.curvesforecasts1[0, 0])
    for path in kept:
        os.remove(path)


@pytest.mark.slow
def test_reproduction_number(fitted):
    import matplotlib.pyplot as plt

    opts, out, _ = fitted
    ropts = sp.options_Rt().replace(mean_GI1=5.0, var_GI1=4.0)
    param_Rt = sp.plotReproductionNumber(AREA, DATE, options=opts, options_Rt=ropts, output_dir=out, verbose=False)
    assert param_Rt.shape == (3, 5) and np.all(param_Rt[:, 0] == 33)  # 30 points and a horizon of 4
    assert np.all(param_Rt[:, 2] <= param_Rt[:, 1]) and np.all(param_Rt[:, 1] <= param_Rt[:, 3])
    assert np.all(param_Rt[:, 1] < 1)  # the curves are past their peak at the end
    assert set(plt.get_fignums()) == {10, 11, 502}
    names_, table = readtable(os.path.join(
        out, "Rt-ranked(1)-onsetfixed-0-flag1-3-method-0-dist-0-horizon-4-daily-coronavirus-deaths-USA-area-52-04-20-2020.csv"))  # fmt: skip
    assert names_ == ["time", "Rt median", "Rt 95%CI LB", "Rt 95% CI UB"]
    assert np.allclose(table[:, 0], np.arange(1, 34)) and np.isnan(table[0, 1]) and np.isfinite(table[-1, 1])
    assert table[-1, 1] == pytest.approx(param_Rt[0, 1])


# ---- the checks of the README -----------------------------------------------------------------


@pytest.mark.slow
def test_quick_start_checks_of_the_readme(tmp_path):
    """The forecast file holds the four later observations, and the baseline has
    a four-day MAE of 63.5; both are independent of the fitted parameters."""
    out = str(tmp_path)
    opts = one_model_options(caddate1="05-11-2020", calibrationperiod1=90)
    sp.Run_Fit_subepidemicFramework(52, "05-11-2020", options=opts, fit_settings=FS, output_dir=out, plot=False,
                                    verbose=False)  # fmt: skip
    ranked, ensembles, baseline, byHorizon = sp.plotForecast_subepidemicFramework(
        52, "05-11-2020", 4, -1, options=opts, output_dir=out, plot=False, verbose=False
    )
    names_, table = readtable(os.path.join(
        out, "ranked(1)-onsetfixed-0-flag1-3-method-0-dist-0-horizon-4-daily-coronavirus-deaths-USA-area-52-05-11-2020.csv"))  # fmt: skip
    assert table[-4:, names_.index("data")].tolist() == [1589, 1735, 1773, 1665]
    assert table.shape[0] == 79 and table[0, :3].tolist() == [2020, 2, 27] and table[-1, :3].tolist() == [2020, 5, 15]
    assert abs(baseline[0] - 63.5) < 1e-8
    assert ensembles.shape == (0, 5) and ranked.shape == (1, 7)
    assert sorted(set(byHorizon[:, 0].tolist())) == [0, 1]


# ---- other data layouts -----------------------------------------------------------------------


@pytest.mark.slow
def test_negative_binomial_fit_writes_alpha(tmp_path):
    out = str(tmp_path)
    opts = one_model_options(method1=3, flag1=sp.GLM, B=5)
    r = sp.Run_Fit_subepidemicFramework(AREA, DATE, options=opts, output_dir=out, plot=False, verbose=False)
    assert opts.dist1 == 3 and r[0].ranked[0].factor1 == pytest.approx(r[0].ranked[0].P0[-2])
    assert np.all(r[0].ranked[0].curves >= 0) and np.all(r[0].ranked[0].curves == np.round(r[0].ranked[0].curves))
    sp.plotFit_subepidemicFramework(AREA, DATE, options=opts, output_dir=out, plot=False, verbose=False)
    names_, table = readtable(files(out, "param-NB-alpha-ranked(1)-*.csv")[0])
    assert names_[0] == "alpha mean" and table.shape == (1, 6) and table[0, 4] == pytest.approx(r[0].ranked[0].P0[-2])
    ranked, _, _, _ = sp.plotForecast_subepidemicFramework(AREA, DATE, 3, -1, options=opts, output_dir=out, plot=False,
                                                           verbose=False)  # fmt: skip
    assert ranked.shape == (1, 7)


@pytest.mark.slow
def test_weekly_cumulative_and_yearly_examples(tmp_path):
    import matplotlib.pyplot as plt

    # weekly, cumulative file, scored against the later file
    out = str(tmp_path / "mpx")
    opts = sp.options_mpx().replace(npatches_fixed=1, numstartpoints=3, B=5)
    r = sp.Run_Fit_subepidemicFramework(options=opts, output_dir=out, plot=False, verbose=False)
    assert r[0].data1.shape[0] == 10
    ranked, _, baseline, byHorizon = sp.plotForecast_subepidemicFramework(
        options=opts, options_forecast=sp.options_forecast().replace(forecastingperiod=3), output_dir=out, verbose=False
    )
    assert ranked.shape == (1, 7) and byHorizon.shape == (6, 7)
    names_, table = readtable(files(out, "ranked(1)-*-horizon-3-weekly-monkeypox-cases-USA-area-1-11-16-2022.csv")[0])
    assert table.shape == (13, 7)
    dn = [sp.datenum(int(y), int(m), int(d)) for y, m, d in table[:, :3]]
    assert set(np.diff(dn).tolist()) == {7}
    plt.close("all")

    # yearly incidence file without later data: not scored
    out = str(tmp_path / "copd")
    opts = sp.options_copd().replace(npatches_fixed=1, numstartpoints=3, B=5)
    sp.Run_Fit_subepidemicFramework(options=opts, output_dir=out, plot=False, verbose=False)
    fopts = sp.options_forecast().replace(getperformance=0, forecastingperiod=3)
    ranked, _, _, _ = sp.plotForecast_subepidemicFramework(options=opts, options_forecast=fopts, output_dir=out, verbose=False)
    assert ranked.shape == (0, 7)
    names_, table = readtable(files(out, "ranked(1)-*-yearly-copd-deaths-CHATHAM-area-1-12-31-2020.csv")[0])
    assert names_ == ["year", "data", "median", "LB", "UB"] and table.shape == (24, 5)
    assert table[:, 0].tolist() == list(range(1, 25)) and np.all(np.isnan(table[21:, 1]))


# ---- standalone curve plots -------------------------------------------------------------------


def test_plot_nsubepidemic_readme_example():
    import matplotlib.pyplot as plt

    totinc = sp.plot_nsubepidemic(1, [0.18, 0.18, 0.18], [0.98, 0.98, 0.98], None, [10000, 5000, 1000], 3, 100, 1, 220, 0)
    assert totinc.shape == (221,) and totinc[0] == pytest.approx(1)
    assert np.sum(totinc) == pytest.approx(16000 - 2, rel=1e-3)  # the three sizes, minus the two seeds
    assert len(plt.gca().lines) == 4
    # the historical rule of the nine-argument call: synchronous when the threshold is zero
    sync = sp.plot_nsubepidemic(1, [0.2, 0.2], [0.9, 0.9], None, [500, 800], 2, 0, 6, 50, plot=False)
    explicit = sp.simulateSubepidemic(np.arange(51.0), 6, 2, 1, 0, 1, [0.2, 0.2], [0.9, 0.9], [1, 1], [500, 800])[2]
    assert np.allclose(sync, explicit)
    # every growth model can be drawn
    for flag1 in range(6):
        curve = sp.plot_nsubepidemic(flag1, [0.3, 0.3], [0.9, 0.9], [1.2, 1.2], [500, 800], 2, 50, 2, 40, 0, plot=False)
        assert curve.shape == (41,) and np.all(np.isfinite(curve))
    with pytest.raises(SubEpiPredictError):
        sp.plot_nsubepidemic(1, None, [0.9], None, [500], 1, 0, 1, 10, plot=False)


def test_plot_nsubepidemic_profiles():
    import matplotlib.pyplot as plt

    curves = sp.plot_nsubepidemic_profiles()
    assert curves.shape == (221, 4) and len(plt.figure(1).axes) == 4
    # a later onset threshold delays the later sub-epidemics
    assert np.argmax(curves[100:, 0]) <= np.argmax(curves[100:, 2])


# ---- rolling forecasts --------------------------------------------------------------------------


def test_rolling_summary_from_hand_made_scores():
    # two dates, baseline (0) and one model (1), horizon 2; columns:
    # date, model, horizon, MAE, MSE, coverage, WIS, WIS of the baseline
    scores = np.array([
        [10, 0, 1, 4.0, 16, 100, 2.0, 2.0], [10, 0, 2, 5.0, 26, 100, 3.0, 3.0],
        [10, 1, 1, 2.0, 4, 100, 1.0, 2.0], [10, 1, 2, 3.0, 10, 50, 2.0, 3.0],
        [17, 0, 1, 8.0, 64, 0, 6.0, 6.0], [17, 0, 2, 8.0, 64, 50, 6.0, 6.0],
        [17, 1, 1, 9.0, 81, 0, 12.0, 6.0], [17, 1, 2, 7.0, 50, 50, 9.0, 6.0],
    ])  # fmt: skip
    summary, columns = sp.subepidemicRollingSummary(scores, 2)
    assert len(columns) == 13 and summary.shape == (4, 13)
    model = summary[(summary[:, 0] == 1) & (summary[:, 1] == 2)][0]
    assert model[2] == 2 and model[3] == 5.0 and model[4] == 50 and model[5] == 5.5
    assert model[6] == pytest.approx(11 / 9)  # ratio of means
    assert model[7] == pytest.approx(np.sqrt((2 / 3) * (9 / 6)))  # geometric mean of the ratios
    assert model[8] == 0.5 and model[9] == pytest.approx(10 / 13)
    # step 2 alone: 2*score(2) - score(1)
    assert model[10] == pytest.approx(((2 * 2 - 1) + (2 * 9 - 12)) / 2)
    assert model[11] == pytest.approx(((2 * 2 - 1) + (2 * 9 - 12)) / ((2 * 3 - 2) + (2 * 6 - 6)))
    assert model[12] == pytest.approx(((2 * 50 - 100) + (2 * 50 - 0)) / 2)
    base = summary[(summary[:, 0] == 0) & (summary[:, 1] == 1)][0]
    assert base[6] == 1 and np.isnan(base[8])
    assert sp.subepidemicRollingSummary(np.zeros((0, 8)), 2)[0].shape == (0, 13)


def test_rolling_forecast_checks_its_inputs(tmp_path):
    kw = dict(options=one_model_options(), output_dir=str(tmp_path), plot=False)
    for args, extra, text in (
        (([],), {}, "Give the forecast dates"),
        ((["2020-04-20"],), {}, "Dates are written mm-dd-yyyy"),
        ((["02-30-2020"],), {}, "not a calendar date"),
        ((["04-20-2020"], 0), {}, "whole number of at least 1"),
        ((["04-20-2020"], 3, [1, 2]), {}, "one area"),
        ((["04-20-2020"], 3, 52, -1, "colour", True), {}, "Unknown option"),
        ((["04-20-2020"], 3, 52, -1, "quiet"), {}, "name and value pairs"),
        ((["01-01-2020"],), {}, "before the first date"),
        ((["04-20-2020"], 3, 99), {}, "52 column(s)"),
        ((["04-20-2020"],), dict(options_forecast=sp.options_forecast().replace(getperformance=0)), "getperformance = 0"),
    ):
        with pytest.raises(SubEpiPredictError) as e:
            sp.subepidemicRollingForecast(*args, **{**kw, **extra})
        assert identifier_of(e) == "SubEpiPredict:RollingInput" and text in e.value.message, text
    # weekly data: a date that is not on the grid of the series
    with pytest.raises(SubEpiPredictError) as e:
        sp.subepidemicRollingForecast(["11-16-2022"], 2, options=sp.options_mpx(), output_dir=str(tmp_path), plot=False)
    assert "not one of the dates of the series" in e.value.message and "11-15-2022 and 11-22-2022" in e.value.message


@pytest.mark.slow
def test_rolling_forecast_scores_dates_and_reuses_fits(tmp_path, capsys):
    out = str(tmp_path / "output")
    inp = str(tmp_path / "input")
    opts = one_model_options()
    kw = dict(options=opts, fit_settings=FS, input_dir=inp, output_dir=out, plot=False)
    # 04-21-2020 has no file of its own: it is written from the most recent file
    dates = ["04-27-2020", "04-20-2020", "04-21-2020"]
    R = sp.subepidemicRollingForecast(dates, 3, AREA, -1, "quiet", True, **kw)
    text = capsys.readouterr().out
    assert R.dates == ["04-20-2020", "04-21-2020", "04-27-2020"] and R.horizon == 3
    assert R.status[:, 1].tolist() == [54, 55, 61]  # data rows on each date
    assert R.status[:, 2].tolist() == [2, 1, 2] and R.status[:, 3].tolist() == [1, 1, 1] and R.status[:, 4].tolist() == [1, 1, 1]
    assert os.path.isfile(os.path.join(inp, "cumulative-daily-coronavirus-deaths-USA-04-21-2020.txt"))
    written = np.loadtxt(os.path.join(inp, "cumulative-daily-coronavirus-deaths-USA-04-21-2020.txt"))
    assert written.shape == (55, 52)
    # scores: 3 dates x (baseline + ranked 1) x 3 horizons
    assert R.scores.shape == (18, 8) and R.summary.shape == (6, 13)
    assert np.all(R.summary[:, 2] == 3) and sorted(set(R.summary[:, 0].tolist())) == [0, 1]
    one = R.scores[(R.scores[:, 1] == 1)]
    assert np.allclose(one[:, 7], R.scores[R.scores[:, 1] == 0][:, 6])  # WIS of the baseline next to every row
    assert os.path.isfile(R.scoresFile) and os.path.isfile(R.summaryFile)
    assert os.path.basename(R.scoresFile) == (
        "rolling-forecast-scores-onsetfixed-0-flag1-3-method-0-dist-0-horizon-3-weight_type--1-daily-coronavirus-deaths"
        "-USA-area-52-04-20-2020-to-04-27-2020.csv")  # fmt: skip
    names_, table = readtable(R.scoresFile)
    assert names_[:5] == ["year", "month", "day", "model", "horizon"] and table.shape == (18, 10)
    assert np.allclose(readtable(R.summaryFile)[1], R.summary, equal_nan=True)
    assert len(files(out, "rolling-log-*.txt")) == 3 and "Rolling forecasts: 3 of 3 date(s) scored" in text
    assert "ratio of means" in text and "Only 3 date(s)" in text

    # a second call uses the saved fits and gives the same scores
    R2 = sp.subepidemicRollingForecast(dates, 3, AREA, -1, quiet=True, **kw)
    assert R2.status[:, 3].tolist() == [2, 2, 2] and np.array_equal(R2.scores, R.scores)
    # other fit settings, or refit=True: the dates are fitted again and the old files are kept aside
    R3 = sp.subepidemicRollingForecast(dates[:1], 3, AREA, -1, quiet=True, **{**kw, "fit_settings": dict(bootstrapStarts=2, **FS)})
    assert R3.status[:, 3].tolist() == [1]
    assert len(os.listdir(os.path.join(out, "replaced-by-rolling"))) == 2  # the screen and rank 1 of 04-27
    R4 = sp.subepidemicRollingForecast(dates[1:2], 3, AREA, -1, refit=True, quiet=True, keepReplaced=False, **kw)
    assert R4.status[:, 3].tolist() == [1]
    # a date whose forecast runs past the data is skipped with a note, not an error
    with pytest.warns(SubEpiPredictWarning, match="RollingDateSkipped"):
        R5 = sp.subepidemicRollingForecast(["05-08-2022"], 3, AREA, -1, quiet=True, **kw)
    assert R5.status[0, 4] == 0 and "skipped" in R5.notes[0] and R5.scoresFile == ""


# ---- every estimation method, error structure and onset mode -----------------------------------

OTHER_SETTINGS = {
    "Poisson likelihood, GLM": dict(method1=1, flag1=sp.GLM, npatches_fixed=1),
    "negative binomial with power d, GLM": dict(method1=5, flag1=sp.GLM, npatches_fixed=1),
    "negative binomial (mean squared), GRM": dict(method1=4, flag1=sp.GRM, npatches_fixed=1),
    "least absolute deviations, Richards": dict(method1=6, dist1=6, flag1=sp.RICH, npatches_fixed=1),
    "least squares with VAR = factor*mean, Gompertz": dict(method1=0, dist1=2, flag1=sp.GOM, npatches_fixed=1),
    "least squares with Poisson noise, GGM": dict(method1=0, dist1=1, flag1=sp.GGM, npatches_fixed=1),
    "synchronous onsets, 2 GLM": dict(method1=0, dist1=0, flag1=sp.GLM, npatches_fixed=2, onset_fixed=1, topmodelsx=2),
    "synchronous onsets, 3 logistic": dict(method1=0, dist1=0, flag1=sp.LM, npatches_fixed=3, onset_fixed=1, topmodelsx=3),
    "negative binomial with a 7-point threshold grid": dict(method1=3, flag1=sp.GLM, npatches_fixed=1, smoothfactor1=7),
}


@pytest.mark.slow
@pytest.mark.parametrize("name", list(OTHER_SETTINGS))
def test_other_settings_run_from_fit_to_reproduction_number(name, tmp_path):
    out = str(tmp_path)
    change = OTHER_SETTINGS[name]
    opts = small_options(**{**dict(numstartpoints=3, B=5), **change})
    kw = dict(options=opts, output_dir=out, plot=False, verbose=False)
    r = sp.Run_Fit_subepidemicFramework(AREA, DATE, fit_settings=FS, **kw)[0]
    best = r.ranked[0]
    n = best.npatches
    assert best.Phatss.shape == (5, 4 * n + 2) and np.all(np.isfinite(best.Phatss)) and np.all(np.isfinite(best.curves))
    method1, dist1 = opts.method1, opts.dist1
    if method1 in (3, 4, 5):
        assert best.factor1 == best.P0[-2] > 0
    elif dist1 in (0, 2, 6):
        assert best.factor1 > 0  # residual scale, or the empirical variance-to-mean ratio
    else:
        assert best.factor1 is None  # Poisson noise has no scale
    if dist1 in (1, 2, 3, 4, 5):
        assert np.all(best.curves >= 0) and np.all(best.curves == np.round(best.curves))
    if opts.onset_fixed == 1:
        assert np.all(r.RMSES[:, 1] == 0) and np.all(np.diff(best.P0[:n]) >= 0)  # one threshold; labels ordered by r
        assert sorted(r.RMSES[:, 0].tolist()) == list(range(1, opts.npatches_fixed + 1))
    performance = sp.plotFit_subepidemicFramework(AREA, DATE, **kw)
    assert performance.shape == (r.topmodelsx, 7) and np.all(np.isfinite(performance))
    assert len(files(out, "param-NB-alpha-*")) == (r.topmodelsx if method1 in (3, 4, 5) else 0)
    assert len(files(out, "param-NB-d-*")) == (r.topmodelsx if method1 == 5 else 0)
    sp.plotRankings_subepidemicFramework(AREA, DATE, **kw)
    ranked, ensembles, baseline, byHorizon = sp.plotForecast_subepidemicFramework(AREA, DATE, 5, 1, **kw)
    assert ranked.shape == (r.topmodelsx, 7) and ensembles.shape == (r.topmodelsx - 1, 5)
    assert np.all(np.isfinite(byHorizon))
    param_Rt = sp.plotReproductionNumber(AREA, DATE, options_Rt=sp.options_Rt().replace(mean_GI1=5.0, var_GI1=4.0), **kw)
    assert param_Rt.shape == (r.topmodelsx, 5) and np.all(np.isfinite(param_Rt)) and np.all(param_Rt[:, 1] > 0)
