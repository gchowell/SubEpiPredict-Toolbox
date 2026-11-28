"""Forecasts of the ranked models and of their ensembles, with scores.

Ports of ``plotForecast_subepidemicFramework.m``,
``getensemblesubepidemicmodels.m``, ``subepidemicAvailableTopModels.m`` and
``subepidemicSkippedRefits.m``.
"""

from __future__ import annotations

import os

import numpy as np

from .doubling import getDoublingTimeCurve
from .errors import SubEpiPredictError, subepidemicWarning, subepidemicWarningLog
from .figures import (
    COMPONENT_COLORS,
    GRAY7,
    GetAdjustedFontSize,
    date_ticks,
    plot_curves,
    series_dates,
    subepidemicPanelGrid,
)
from .io import QUANT_NAMES_RANKED, OutputNames, getData, load_mat, save_mat, scalar, vector, writetable
from .model import subepidemicSimulateRefit
from .noise import AddPoissonError
from .optionsets import options as default_options
from .optionsets import options_forecast as default_options_forecast
from .optionsets import resolve_options
from .scores import (
    _clip0,
    computeforecastperformance,
    computeQuantiles,
    computeWIS,
    subepidemicBaselineForecast,
    subepidemicEnsembleWeights,
    subepidemicWeightedQuantile,
)
from .settings import resolve_fit_settings, subepidemicSeed
from .search import unique_rows_stable as unique_rows
from .utils import cadtemporal_of, datenum, datevec, num2ordinal, num2str, quantile

DOUBLING_NAMES = ["i_th doubling", "db mean", "db 95%CI LB", "db 95% CI UB", "prob. i_th doubling"]
FORECAST_NAMES = ["year", "month", "day", "data", "median", "LB", "UB"]


# <============================================================================>
# <================================ Small helpers =============================>
# <============================================================================>


def subepidemicAvailableTopModels(topmodels, names):
    """Limit the ranked models to those the saved fit holds.

    Port of ``subepidemicAvailableTopModels.m``. ``topmodels`` holds the ranks
    a function is about to process (1 to ``topmodelsx``); ``names`` is the
    :class:`~subepipredict.io.OutputNames` of the area and date. Returns the
    leading part for which a ranked file of the fit exists, with warning
    ``SubEpiPredict:TopModelsCapped`` when it is shorter. If no ranked file
    exists at all, the list is returned unchanged and the caller stops at its
    own load.
    """
    topmodels = [int(k) for k in topmodels]
    nAvailable = 0
    for k in topmodels:
        if os.path.isfile(names.ranked_mat(k)):
            nAvailable += 1
        else:
            break
    if 0 < nAvailable < len(topmodels):
        subepidemicWarning(
            "SubEpiPredict:TopModelsCapped",
            "The options ask for %d top-ranked models, but the saved fit has %d. Ranks 1 to %d are used.",
            len(topmodels), nAvailable, nAvailable,
        )  # fmt: skip
        topmodels = topmodels[:nAvailable]
    return topmodels


def subepidemicSkippedRefits(label, rank1, nSkipped, nRefits, firstMessage, settings=None):
    """One policy for refits that could not be simulated.

    Port of ``subepidemicSkippedRefits.m``. ``nSkipped`` of the ``nRefits``
    could not be simulated over the forecast period and were left out. None
    skipped: nothing happens. All skipped: error
    ``SubEpiPredict:NoForecastCurves``. Share above the setting
    ``maxSkippedRefitFraction`` (default 0.05): error
    ``SubEpiPredict:TooManyRefitsSkipped``. Otherwise: warning
    ``SubEpiPredict:ForecastCurvesSkipped``.
    """
    if nSkipped <= 0:
        return
    settings = resolve_fit_settings(settings)
    limit = settings.maxSkippedRefitFraction
    nUsed = nRefits - nSkipped
    if nUsed <= 0:
        raise SubEpiPredictError(
            "SubEpiPredict:NoForecastCurves",
            "%s, ranked model %d: none of the %d refits could be simulated over the forecast period. "
            "First failure: %s" % (label, rank1, nRefits, firstMessage),
        )
    if nSkipped / nRefits > limit:
        raise SubEpiPredictError(
            "SubEpiPredict:TooManyRefitsSkipped",
            "%s, ranked model %d: %d of %d refits could not be simulated over the forecast period. "
            "That is more than the share allowed (maxSkippedRefitFraction = %g). First failure: %s"
            % (label, rank1, nSkipped, nRefits, limit, firstMessage),
        )
    subepidemicWarning(
        "SubEpiPredict:ForecastCurvesSkipped",
        "%s, ranked model %d: %d of %d refits could not be simulated over the forecast period and were left out. "
        "The result rests on %d refits. First failure: %s",
        label, rank1, nSkipped, nRefits, nUsed, firstMessage,
    )  # fmt: skip


def load_ranked(names, rank):
    """The saved fit of one ranked model as a Struct of plain values.

    Reads ``modifiedLogisticPatch-ensem-...-rank-<k>.mat`` and unwraps the
    variables the fit-plot, forecast and Rt functions use.
    """
    path = names.ranked_mat(rank)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            "The saved fit %s was not found. Run Run_Fit_subepidemicFramework with the same options first." % path
        )
    raw = load_mat(path)
    npatches = int(scalar(raw.npatches))
    nz = 4 * npatches + 2
    S = raw.copy()
    S.npatches = npatches
    S.onset_thr = scalar(raw.onset_thr)
    S.AICc_best = scalar(raw.AICc_best)
    S.data1 = np.asarray(raw.data1, dtype=float).reshape(-1, 2)
    S.data = np.asarray(raw.data, dtype=float).reshape(-1, 2) if "data" in raw else S.data1
    S.Phatss = np.asarray(raw.Phatss, dtype=float).reshape(-1, nz)
    S.curves = np.asarray(raw.curves, dtype=float)
    if S.curves.shape[0] != S.data1.shape[0]:
        S.curves = S.curves.reshape(S.data1.shape[0], -1)
    S.bestfit = vector(raw.bestfit) if "bestfit" in raw else None
    S.P0 = vector(raw.P0) if "P0" in raw else None
    S.relativelik_i = vector(raw.relativelik_i)
    S.t_window = vector(raw.t_window).astype(int)
    S.timelags = int(scalar(raw.timelags))
    S.flag1 = int(vector(raw.flag1)[0])
    S.M = S.Phatss.shape[0]
    f = np.asarray(raw.get("factor1", np.zeros(0)), dtype=float).ravel()
    S.factor1 = float(f[0]) if f.size else None
    S.rank1 = int(scalar(raw.rank1)) if "rank1" in raw else int(rank)
    S.outbreak1 = int(scalar(raw.outbreak1)) if "outbreak1" in raw else int(names.outbreakx)
    return S


def _doubling_summary(curves1, DT, M, weights=None):
    """Sequence of doubling times over the curves (columns) and their mean.

    Returns ``(seq_doublingtimes, param_doubling)``: rows ``[i, mean, 95% LB,
    95% UB, share of curves with an i-th doubling]`` and ``[mean, LB, UB]`` of
    the mean doubling time of a curve. With ``weights`` every curve enters
    with its weight (ensembles); the share is then a share of the weight.
    """
    n = curves1.shape[1]
    meandoublingtime = np.full(n, np.nan)
    doublingtimess = np.full((30, n), np.nan)
    maxd = 1
    for j in range(n):
        _, _, _, doublingtimes = getDoublingTimeCurve(np.maximum(curves1[:, j], 0), DT, 0)
        doublingtimess[: doublingtimes.size, j] = doublingtimes
        maxd = max(maxd, doublingtimes.size)
        if doublingtimes.size:
            meandoublingtime[j] = np.mean(doublingtimes)
    doublingtimess = doublingtimess[:maxd, :]

    seq = []
    for j in range(maxd):
        index1 = np.flatnonzero(~np.isnan(doublingtimess[j, :]))
        if weights is None:
            if index1.size == 0:
                seq.append([j + 1, np.nan, np.nan, np.nan, 0.0])
            else:
                v = doublingtimess[j, index1]
                seq.append([j + 1, np.mean(v), quantile(v, 0.025), quantile(v, 0.975), index1.size / M])
        else:
            if index1.size == 0:
                seq.append([j + 1, np.nan, np.nan, np.nan, 0.0])
            else:
                wj = weights[index1]
                v = doublingtimess[j, index1]
                lim = subepidemicWeightedQuantile(v, [0.025, 0.975], wj)[0]
                seq.append([j + 1, np.sum(wj * v) / np.sum(wj), lim[0], lim[1], np.sum(wj) / np.sum(weights)])
    seq = np.array(seq, dtype=float).reshape(-1, 5)

    index1 = np.flatnonzero(~np.isnan(meandoublingtime))
    if index1.size == 0:  # no curve doubled
        param_doubling = np.array([np.nan, np.nan, np.nan])
    elif weights is None:
        v = meandoublingtime[index1]
        param_doubling = np.array([np.mean(v), quantile(v, 0.025), quantile(v, 0.975)])
    else:
        wj = weights[index1]
        v = meandoublingtime[index1]
        lim = subepidemicWeightedQuantile(v, [0.025, 0.975], wj)[0]
        param_doubling = np.array([np.sum(wj * v) / np.sum(wj), lim[0], lim[1]])
    return seq, param_doubling


def _forecast_table(datevecfirst1, timelags, DT, observed, median, LB1, UB1):
    """Rows ``[year, month, day, data, median, LB, UB]`` of a forecast file."""
    n = len(median)
    dn = series_dates(datevecfirst1, timelags, DT, n)
    ymd = np.array([datevec(v) for v in dn], dtype=float)
    return np.column_stack([ymd, observed, median, LB1, UB1])


def _interval_panel(ax, timevect2, LB1, UB1, median, data1, timevect_end, fs):
    """Gray 95% band, median, limits, data and the line at the end of the calibration period."""
    ax.fill_between(timevect2, LB1, UB1, color=(0.8, 0.8, 0.8))
    ax.plot(timevect2, median, "r-", linewidth=2)
    ax.plot(timevect2, LB1, "k--", linewidth=2)
    ax.plot(timevect2, UB1, "k--", linewidth=2)
    ax.plot(data1[:, 0], data1[:, 1], "ko", markerfacecolor="none", markeredgewidth=1.5)
    top = float(np.nanmax(UB1)) * 1.2 if np.nanmax(UB1) > 0 else 1.0
    ax.set_xlim(0, len(timevect2) - 1)
    ax.set_ylim(0, top)
    ax.plot([timevect_end, timevect_end], [0, top], "k--", linewidth=2)
    ax.tick_params(labelsize=fs - 3)


# <============================================================================>
# <================================= Ensembles ================================>
# <============================================================================>


def getensemblesubepidemicmodels(
    cumulative1, cadfilename1, datevecfirst1, npatches_fixed, onset_fixed, smoothfactor1, outbreakx, cadregion,
    caddate1, caddisease, datatype, flag1, method1, dist1, calibrationperiod1, topmodels1, forecastingperiod,
    getperformance, weight_type1, WISC_hash, WISF_hash, printscreen1, ensemblePanelGrid=None, *,
    input_dir="./input", output_dir="./output",
):  # fmt: skip
    """Ensemble forecast of the top-ranked models ``topmodels1`` and its scores.

    Port of ``getensemblesubepidemicmodels.m``. The ensemble is the mixture of
    its members: member k with weight ``w(k)``
    (:func:`subepidemicEnsembleWeights`) gives ALL its curves, each with
    weight ``w(k)/M``; medians, intervals, quantiles and scores of the
    ensemble are weighted quantiles (:func:`subepidemicWeightedQuantile`).

    The function reads the ranking (``ABC-ensem-...mat``) and the forecast
    curves of the members (``Forecast-modifiedLogisticPatch-...-rank-<k>.mat``,
    written by :func:`plotForecast_subepidemicFramework`) from ``output_dir``
    and writes ``quantileTimes-Ensemble(k)-``, ``doublingTimes-Ensemble(k)-``
    and ``Ensemble(k)-...csv``.

    The arguments are those of the MATLAB function; ``topmodels1`` holds the
    ONE-BASED ranks of the members (1 to k), ``WISC_hash`` and ``WISF_hash``
    the scores of the ranked models for the weight codes 2 and 3.

    Returns ``(RMSECS, MSECS, MAECS, PICS, MISCS, WISC, RMSEFS, MSEFS, MAEFS,
    PIFS, MISFS, WISFS, forecast1, quantilesc, quantilesf)``; ``forecast1``
    has the columns median, upper and lower 95% limit.
    """
    topmodels1 = [int(k) for k in np.asarray(topmodels1).ravel()]
    last = topmodels1[-1]

    abcName = (
        "ABC-ensem-npatchesfixed-%s-onsetfixed-%s-smoothing-%s-%s-flag1-%s-method-%s-dist-%s-calibrationperiod-%s.mat"
        % (num2str(npatches_fixed), num2str(onset_fixed), num2str(smoothfactor1), cadfilename1, num2str(flag1),
           num2str(method1), num2str(dist1), num2str(calibrationperiod1))
    )  # fmt: skip
    abc = load_mat(os.path.join(output_dir, abcName))
    RMSES = np.asarray(abc.RMSES, dtype=float).reshape(-1, 3)
    data1 = np.asarray(abc.data1, dtype=float).reshape(-1, 2)
    timevect = data1[:, 0]
    DT = int(scalar(abc.DT))
    cadtemporal = cadtemporal_of(DT)

    # remove any repeated rows: [npatches onset_thr AICc]
    RMSES, _ = unique_rows(RMSES)

    # Weights of the member models: -1 equal, 0 and 1 Akaike weights, 2 inverse
    # WISC, 3 inverse WISF. The scores are read only for the code that uses them.
    scoresC = None
    scoresF = None
    index = np.array(topmodels1) - 1
    if weight_type1 == 2:
        scoresC = np.asarray(WISC_hash, dtype=float).ravel()[index]
    elif weight_type1 == 3:
        scoresF = np.asarray(WISF_hash, dtype=float).ravel()[index]
    weights1 = subepidemicEnsembleWeights(weight_type1, RMSES[index, 2], scoresC, scoresF)

    ens1, ens2, w1, w2 = [], [], [], []
    timevect2 = None
    timelags = 0
    for position, rank1 in enumerate(topmodels1):
        forecastName = (
            "Forecast-modifiedLogisticPatch-ensem-npatchesfixed-%s-onsetfixed-%s-smoothing-%s-%s-flag1-%s-method-%s"
            "-dist-%s-calibrationperiod-%s-forecastingperiod-%s-rank-%s.mat"
            % (num2str(npatches_fixed), num2str(onset_fixed), num2str(smoothfactor1), cadfilename1, num2str(flag1),
               num2str(method1), num2str(dist1), num2str(calibrationperiod1), num2str(forecastingperiod),
               num2str(rank1))
        )  # fmt: skip
        F = load_mat(os.path.join(output_dir, forecastName))
        curvesforecasts1 = np.asarray(F.curvesforecasts1, dtype=float)
        curvesforecasts2 = np.asarray(F.curvesforecasts2, dtype=float)
        timevect2 = vector(F.timevect2)
        timelags = int(scalar(F.timelags))
        datevecend1 = [int(v) for v in vector(F.datevecend1)]
        M1 = curvesforecasts1.shape[1]
        M2 = curvesforecasts2.shape[1]
        if weights1[position] > 0 and M1 > 0 and M2 > 0:
            ens1.append(curvesforecasts1)
            w1.append(np.full(M1, weights1[position] / M1))
            ens2.append(curvesforecasts2)
            w2.append(np.full(M2, weights1[position] / M2))

    if not ens1 or not ens2:
        raise SubEpiPredictError(
            "SubEpiPredict:EnsembleEmpty",
            "Ensemble of %d models: no member with a positive weight has forecast curves." % len(topmodels1),
        )
    curvesforecasts1 = np.hstack(ens1)
    curvesforecasts2 = np.hstack(ens2)
    # the weights sum to 1 (they already do unless a member was left out above)
    ensWeights1 = np.concatenate(w1)
    ensWeights1 = ensWeights1 / ensWeights1.sum()
    ensWeights2 = np.concatenate(w2)
    ensWeights2 = ensWeights2 / ensWeights2.sum()

    # Median and 95% limits of the ensemble.
    ensLimits = subepidemicWeightedQuantile(curvesforecasts2, [0.025, 0.5, 0.975], ensWeights2)
    ensMedian = ensLimits[:, 1]
    LB1 = _clip0(ensLimits[:, 0])
    UB1 = _clip0(ensLimits[:, 2])

    # store forecast curves
    forecast1 = np.column_stack([ensMedian, ensLimits[:, 2], ensLimits[:, 0]])

    quantilesc, quantilesf, _ = computeQuantiles(data1, curvesforecasts2, forecastingperiod, ensWeights2)

    names = _EnsembleNames(output_dir, onset_fixed, flag1, method1, dist1, cadtemporal, caddisease, datatype,
                           cadregion, outbreakx, caddate1)  # fmt: skip

    # Exporting the quantile forecasts
    writetable(
        names.csv("quantileTimes-Ensemble(%s)" % num2str(last), forecastingperiod, weight_type1, "weighttype"),
        np.vstack([quantilesc, quantilesf]), QUANT_NAMES_RANKED,
    )  # fmt: skip

    ax = None
    if printscreen1:
        import matplotlib.pyplot as plt

        fs = GetAdjustedFontSize()
        # One panel per ensemble, Ensemble(2) in panel 1.
        if ensemblePanelGrid is None:
            ensemblePanelGrid = [1, max(3, last - 1)]
        if ensemblePanelGrid[0] * ensemblePanelGrid[1] < last - 1:
            ensemblePanelGrid = [int(np.ceil((last - 1) / 3)), 3]
        if plt.fignum_exists(500):
            fig = plt.figure(500)
        else:
            fig = plt.figure(500, figsize=(6.2 * ensemblePanelGrid[1], 4.6 * ensemblePanelGrid[0]))
        ax = fig.add_subplot(ensemblePanelGrid[0], ensemblePanelGrid[1], last - 1)
        _interval_panel(ax, timevect2, LB1, UB1, ensMedian, data1, timevect[-1], fs)
        date_ticks(ax, datevecfirst1, timelags, DT, len(timevect2))
        ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
        ax.set_title("Ensemble(%s)" % num2str(last), fontsize=fs, fontweight="bold")
        fig.set_facecolor("white")

    # compute doubling times: every curve of every member enters with its weight
    seq_doublingtimes, _ = _doubling_summary(curvesforecasts1, DT, curvesforecasts1.shape[1], ensWeights1)

    # <=============================================================================================>
    # <============================== Save file with doubling time estimates =======================>
    # <=============================================================================================>
    writetable(
        names.csv("doublingTimes-Ensemble(%s)" % num2str(last), forecastingperiod, weight_type1, "weighttype"),
        seq_doublingtimes, DOUBLING_NAMES,
    )  # fmt: skip

    forecastFile = names.csv("Ensemble(%s)" % num2str(last), forecastingperiod, weight_type1, "weighttype")
    nan = np.nan
    empty = np.zeros((0, 2))
    if getperformance:
        # most recent data
        datenum1 = datenum(caddate1)
        if DT != 365:
            datenum1 = datenum1 + DT
        horizonScored = int(forecastingperiod * getperformance)
        if horizonScored > 0:
            data2 = getData(cumulative1, cadtemporal, caddisease, datatype, cadregion, DT, datevecfirst1,
                            datevecend1, datevec(datenum1), outbreakx, horizonScored, input_dir)  # fmt: skip
            times2 = data1[-1, 0] + 1 + np.arange(horizonScored)
            if ax is not None:
                ax.plot(times2, data2, "ro", markerfacecolor="none", markeredgewidth=1.5)
            datalatest2 = np.vstack([data1, np.column_stack([times2, data2])])
            observed = np.concatenate([data1[:, 1], data2])
        else:
            datalatest2 = data1
            observed = data1[:, 1]
        # <=============================================================================================>
        # <============================== Save file with forecast ======================================>
        # <=============================================================================================>
        writetable(forecastFile, _forecast_table(datevecfirst1, timelags, DT, observed, ensMedian, LB1, UB1),
                   FORECAST_NAMES)  # fmt: skip

        # compute performance metrics
        (RMSECS, MSECS, MAECS, PICS, MISCS, RMSEFS, MSEFS, MAEFS, PIFS, MISFS) = computeforecastperformance(
            data1, datalatest2, curvesforecasts1, curvesforecasts2, horizonScored, ensWeights1, ensWeights2
        )
        WISC, WISFS, _, _ = computeWIS(data1, datalatest2, curvesforecasts2, horizonScored, ensWeights2)
    else:
        RMSECS = MSECS = MAECS = PICS = MISCS = WISC = nan
        RMSEFS = MSEFS = MAEFS = PIFS = MISFS = WISFS = empty
        observed = np.concatenate([data1[:, 1], np.full(int(forecastingperiod), np.nan)])
        writetable(forecastFile, _forecast_table(datevecfirst1, timelags, DT, observed, ensMedian, LB1, UB1),
                   FORECAST_NAMES)  # fmt: skip

    return (RMSECS, MSECS, MAECS, PICS, MISCS, WISC, RMSEFS, MSEFS, MAEFS, PIFS, MISFS, WISFS, forecast1,
            quantilesc, quantilesf)  # fmt: skip


class _EnsembleNames:
    """CSV names of :func:`getensemblesubepidemicmodels` from its own arguments."""

    def __init__(self, output_dir, onset_fixed, flag1, method1, dist1, cadtemporal, caddisease, datatype, cadregion,
                 outbreakx, caddate1):  # fmt: skip
        self.output_dir = output_dir
        self.onset_fixed = onset_fixed
        self.flag1 = flag1
        self.method1 = method1
        self.dist1 = dist1
        self.cadtemporal = cadtemporal
        self.caddisease = caddisease
        self.datatype = datatype
        self.cadregion = cadregion
        self.outbreakx = outbreakx
        self.caddate1 = caddate1

    csv = OutputNames.csv
    _path = OutputNames._path


# <============================================================================>
# <============================ The forecast function =========================>
# <============================================================================>


def plotForecast_subepidemicFramework(
    outbreakx_pass=None,
    caddate1_pass=None,
    forecastingperiod_pass=None,
    weight_type1_pass=None,
    options=None,
    options_forecast=None,
    fit_settings=None,
    input_dir="./input",
    output_dir="./output",
    plot=True,
    verbose=True,
):
    """Forecast with the ranked models and their ensembles, and score the forecasts.

    Port of ``plotForecast_subepidemicFramework.m``. Run
    :func:`Run_Fit_subepidemicFramework` with the same options first.

    Parameters
    ----------
    outbreakx_pass, caddate1_pass : area and date stamp (None = the options)
    forecastingperiod_pass : forecast horizon (None = ``options_forecast``)
    weight_type1_pass : ensemble weights: -1 equal, 1 Akaike, 2 inverse
        calibration WIS, 0 as 1 with a warning, 3 as 2 (None = ``options_forecast``)
    options, options_forecast : option objects or functions returning them
    fit_settings : fit settings (random seed, skipped-refit limit, baseline window)
    input_dir, output_dir : folders of the data and of the results
    plot : draw the figures (10 ranked forecasts with components, 11 ranked
        forecasts with intervals, 500 ensembles, 400 and 401 scores)

    Every bootstrap parameter vector of a ranked model gives one model
    trajectory over the calibration and forecast period and 20 simulated
    observation series around it. Ensemble(k) combines ranked models 1 to k.
    With ``getperformance = 1`` the forecasts are scored against the later
    observations of the data file dated ``datevecend1``, next to a flat
    baseline (:func:`subepidemicBaselineForecast`).

    Returns
    -------
    ``(performanceTop, performanceEns, performanceBaseline, scoresByHorizonOut)``
        ``performanceTop``: one row per ranked model, ``[rank, MAE, MSE,
        coverage of the 95% PI, WIS, AICc, relative likelihood]``.
        ``performanceEns``: one row per ensemble, ``[size, MAE, MSE, coverage,
        WIS]``. ``performanceBaseline``: ``[MAE, MSE, coverage, WIS]`` of the
        baseline over the full horizon. ``scoresByHorizonOut``: one row per
        model and horizon, ``[model, horizon, MAE, MSE, coverage, WIS, WIS
        relative to the baseline]``; model 0 is the baseline, k ranked model
        k, 100+k Ensemble(k); the scores of horizon h average the first h
        time units. All are empty arrays when ``getperformance = 0``.

    Files (names as in MATLAB): ``quantile-ranked(k)-``,
    ``doublingTimes-ranked(k)-``, ``ranked(k)-``, ``quantileTimes-Ensemble(k)-``,
    ``doublingTimes-Ensemble(k)-``, ``Ensemble(k)-``,
    ``performance-forecasting-topRanked-``, ``performance-forecasting-Ensemble-``,
    ``performance-forecasting-byHorizon-``, ``performance-calibration-topRanked-``
    and ``performance-calibration-Ensemble-...csv``.
    """
    opts = resolve_options(options, default_options)
    fopts = resolve_options(options_forecast, default_options_forecast)
    fitSettings = resolve_fit_settings(fit_settings)

    performanceTop = np.zeros((0, 7))
    performanceEns = np.zeros((0, 5))
    performanceBaseline = np.zeros(0)
    scoresByHorizonOut = np.zeros((0, 7))

    outbreakx = opts.outbreakx if outbreakx_pass is None or np.size(outbreakx_pass) == 0 else outbreakx_pass
    outbreakx = int(np.asarray(outbreakx).ravel()[0])
    caddate1 = opts.caddate1 if caddate1_pass is None or caddate1_pass == "" else caddate1_pass

    cumulative1 = opts.cumulative1
    cadregion = opts.cadregion
    caddisease = opts.caddisease
    datatype = opts.datatype
    datevecfirst1 = opts.datevecfirst1
    datevecend1 = opts.datevecend1
    DT = opts.DT
    cadtemporal = cadtemporal_of(DT)
    npatches_fixed = opts.npatches_fixed
    onset_fixed = opts.onset_fixed
    method1 = opts.method1
    dist1 = opts.dist1
    smoothfactor1 = opts.smoothfactor1
    calibrationperiod1 = opts.calibrationperiod1

    topmodels1 = list(range(1, opts.topmodelsx + 1))
    if npatches_fixed == 1:
        topmodels1 = [1]

    names = OutputNames(opts, outbreakx, caddate1, output_dir)
    cadfilename2 = names.cadfilename2

    # A fresh list of warnings for this run.
    subepidemicWarningLog("start")

    # The options can ask for more top models than the fit saved; the ranks
    # that exist are processed.
    topmodels1 = subepidemicAvailableTopModels(topmodels1, names)

    # Random seed (setting rngSeed): the forecast noise repeats from run to run.
    subepidemicSeed("plotForecast_subepidemicFramework", fitSettings, verbose=verbose)

    rows, cols = subepidemicPanelGrid(len(topmodels1))
    # Grid of the ensemble figure: one panel per ensemble, three per row.
    ensemblePanelGrid = [max(1, int(np.ceil((len(topmodels1) - 1) / 3))), 3]

    getperformance = int(fopts.getperformance)
    deletetempfiles = int(fopts.deletetempfiles)
    forecastingperiod = fopts.forecastingperiod if forecastingperiod_pass is None else forecastingperiod_pass
    forecastingperiod = int(forecastingperiod)
    printscreen1 = bool(plot)
    weight_type1 = fopts.weight_type1 if weight_type1_pass is None else weight_type1_pass

    nTop = max(topmodels1)
    WISC_hash = np.zeros(nTop)  # WISC of the calibration period, for weight_type1 = 2
    WISF_hash = np.zeros(nTop)  # WISF of the forecast, kept as in MATLAB

    # rows [model, area, date number, value]
    RMSECSS, MSECSS, MAECSS, PICSS, MISCSS, WISCSS = [], [], [], [], [], []
    RMSEFSS, MSEFSS, MAEFSS, PIFSS, MISFSS, WISFSS = [], [], [], [], [], []
    scoresByHorizon = []
    AICc_rank1 = []
    relativelik_rank1 = []
    datalatest2 = None
    data1 = None
    scored = getperformance == 1 and forecastingperiod > 0
    dnum = datenum(caddate1)

    if printscreen1:
        import matplotlib.pyplot as plt

        plt.close("all")
        fs = GetAdjustedFontSize()
        fig10 = plt.figure(10, figsize=(6.4 * cols, 4.4 * rows))
        fig11 = plt.figure(11, figsize=(6.4 * cols, 4.4 * rows))

    for cc1, rankx in enumerate(topmodels1, start=1):

        S = load_ranked(names, rankx)
        npatches = S.npatches
        onset_thr = S.onset_thr
        flag1 = S.flag1
        Phatss = S.Phatss
        data1 = S.data1
        timelags = S.timelags
        rank1 = S.rank1
        outbreak1 = S.outbreak1
        factor1 = S.factor1
        M = S.M

        AICc_rank1.append([rank1, S.AICc_best])
        relativelik_rank1.append([rank1, S.relativelik_i[rankx - 1]])

        d_hat = 1

        timevect = data1[:, 0]
        timevect2 = np.arange(0, S.t_window[-1] + forecastingperiod, dtype=float)

        cols1 = []
        cols2 = []
        components = [[] for _ in range(npatches)]

        nSimulationsSkipped = 0  # refits left out because their curve could not be simulated
        firstSkipMessage = ""

        # <======================================================================================>
        # <== generate forecast curves from each bootstrap realization ==========================>
        # <======================================================================================>
        for realization in range(M):
            rs_hat = Phatss[realization, 0:npatches]
            ps_hat = Phatss[realization, npatches : 2 * npatches]
            as_hat = Phatss[realization, 2 * npatches : 3 * npatches]
            Ks_hat = Phatss[realization, 3 * npatches : 4 * npatches]

            if method1 == 3 or method1 == 4:
                factor1 = Phatss[realization, -2]
            elif method1 == 5:
                factor1 = Phatss[realization, -2]
                d_hat = Phatss[realization, -1]

            # A refit whose curve cannot be integrated over the forecast period
            # is left out and counted.
            simulated, x, totinc, patchIncidence, skipMessage = subepidemicSimulateRefit(
                timevect2, data1[0, 1], npatches, onset_fixed, onset_thr, flag1, rs_hat, ps_hat, as_hat, Ks_hat
            )
            if not simulated:
                nSimulationsSkipped += 1
                if not firstSkipMessage:
                    firstSkipMessage = "refit %d: %s" % (realization + 1, skipMessage)
                continue

            for j in range(npatches):
                components[j].append(patchIncidence[:, j])

            cols1.append(totinc)
            forecasts2 = AddPoissonError(np.cumsum(totinc), 20, dist1, factor1, d_hat)
            cols2.append(forecasts2)

        nForecastCurves = len(cols1)

        # none skipped: nothing; all skipped or too many: error; otherwise a warning
        subepidemicSkippedRefits("Forecast", rank1, nSimulationsSkipped, M, firstSkipMessage, fitSettings)

        curvesforecasts1 = np.column_stack(cols1)
        curvesforecasts2 = np.hstack(cols2)

        quantilesc, quantilesf, _ = computeQuantiles(data1, curvesforecasts2, forecastingperiod)

        # Exporting the quantile forecasts
        writetable(names.csv("quantile-ranked(%s)" % num2str(rank1), forecastingperiod),
                   np.vstack([quantilesc, quantilesf]), QUANT_NAMES_RANKED)  # fmt: skip

        # compute doubling times
        seq_doublingtimes, param_doubling = _doubling_summary(curvesforecasts1, DT, M)

        # <=============================================================================================>
        # <============================== Save file with doubling time estimates =======================>
        # <=============================================================================================>
        writetable(names.csv("doublingTimes-ranked(%s)" % num2str(rank1), forecastingperiod),
                   seq_doublingtimes, DOUBLING_NAMES)  # fmt: skip

        if printscreen1:
            ax = fig10.add_subplot(rows, cols, cc1)
            for j in range(npatches):
                plot_curves(ax, timevect2, np.column_stack(components[j]), COMPONENT_COLORS[j % 7])
            plot_curves(ax, timevect2, curvesforecasts1, GRAY7)
            ax.set_title("%s Ranked Model" % num2ordinal(rank1), fontsize=fs, fontweight="bold")
            ax.plot(data1[:, 0], data1[:, 1], "ko", markerfacecolor="none", markeredgewidth=1.5)
            top = float(np.max(data1[:, 1])) * 2
            ax.set_xlim(0, len(timevect2) - 1)
            ax.set_ylim(0, top)
            ax.plot([timevect[-1], timevect[-1]], [0, top], "k--", linewidth=2)
            date_ticks(ax, datevecfirst1, timelags, DT, len(timevect2))
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.tick_params(axis="y", labelsize=fs - 3)
            fig10.set_facecolor("white")

        median2 = np.median(curvesforecasts2, axis=1)
        LB1 = quantile(curvesforecasts2, 0.025, axis=1)
        LB1 = _clip0(LB1)
        UB1 = quantile(curvesforecasts2, 0.975, axis=1)
        UB1 = _clip0(UB1)

        ax11 = None
        if printscreen1:
            ax11 = fig11.add_subplot(rows, cols, cc1)
            _interval_panel(ax11, timevect2, LB1, UB1, median2, data1, timevect[-1], fs)
            date_ticks(ax11, datevecfirst1, timelags, DT, len(timevect2))
            ax11.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax11.set_title("%s Ranked Model" % num2ordinal(rank1), fontsize=fs, fontweight="bold")
            fig11.set_facecolor("white")

        # save forecast realizations of the model
        save_mat(
            names.forecast_mat(forecastingperiod, rankx),
            dict(curvesforecasts1=curvesforecasts1, curvesforecasts2=curvesforecasts2, datevecfirst1=datevecfirst1,
                 datevecend1=datevecend1, timevect2=timevect2, timelags=timelags, cadtemporal=cadtemporal,
                 nSimulationsSkipped=nSimulationsSkipped, nForecastCurves=nForecastCurves),
        )  # fmt: skip

        # <=========================================================================================>
        # <================================ Load forecast period data ==============================>
        # <=========================================================================================>
        datenum1 = dnum
        if DT != 365:
            datenum1 = datenum1 + DT

        nRows = data1.shape[0] + forecastingperiod
        if scored:
            data2 = getData(cumulative1, cadtemporal, caddisease, datatype, cadregion, DT, datevecfirst1,
                            datevecend1, datevec(datenum1), outbreak1, forecastingperiod, input_dir)  # fmt: skip
            times2 = data1[-1, 0] + 1 + np.arange(forecastingperiod)
            if ax11 is not None:
                ax11.plot(times2, data2, "ro", markerfacecolor="none", markeredgewidth=1.5)
            observed = np.concatenate([data1[:, 1], data2])
            datalatest2 = np.vstack([data1, np.column_stack([times2, data2])])
            rankedFile = names.csv("ranked(%s)" % num2str(rank1), forecastingperiod)
        else:
            observed = np.concatenate([data1[:, 1], np.full(forecastingperiod, np.nan)])
            datalatest2 = data1
            rankedFile = names.csv("ranked(%s)" % num2str(rank1))

        # <=============================================================================================>
        # <============================== Save file with forecast ======================================>
        # <=============================================================================================>
        if DT == 365:
            writetable(rankedFile,
                       np.column_stack([np.arange(1, nRows + 1), observed, median2, LB1, UB1]),
                       ["year", "data", "median", "LB", "UB"])  # fmt: skip
        else:
            writetable(rankedFile, _forecast_table(datevecfirst1, timelags, DT, observed, median2, LB1, UB1),
                       FORECAST_NAMES)  # fmt: skip

        # <=============================================================================================>
        # <================================== Performance metrics ======================================>
        # <=============================================================================================>
        horizonScored = forecastingperiod * getperformance
        (RMSECS_model1, MSECS_model1, MAECS_model1, PICS_model1, MISCS_model1, RMSEFS_model1, MSEFS_model1,
         MAEFS_model1, PIFS_model1, MISFS_model1) = computeforecastperformance(
            data1, datalatest2, curvesforecasts1, curvesforecasts2, horizonScored
        )  # fmt: skip
        WISC, WISFS, _, _ = computeWIS(data1, datalatest2, curvesforecasts2, horizonScored)

        # WISC of the model, for the weights of the ensemble (weight_type1 = 2)
        WISC_hash[rankx - 1] = WISC
        if scored:
            WISF_hash[rankx - 1] = WISFS[-1, -1]

        key = [rankx, outbreakx, dnum]
        RMSECSS.append(key + [RMSECS_model1])
        MSECSS.append(key + [MSECS_model1])
        MAECSS.append(key + [MAECS_model1])
        PICSS.append(key + [PICS_model1])
        MISCSS.append(key + [MISCS_model1])
        WISCSS.append(key + [WISC])

        if scored:
            RMSEFSS.append(key + [RMSEFS_model1[-1, -1]])
            MSEFSS.append(key + [MSEFS_model1[-1, -1]])
            MAEFSS.append(key + [MAEFS_model1[-1, -1]])
            PIFSS.append(key + [PIFS_model1[-1, -1]])
            MISFSS.append(key + [MISFS_model1[-1, -1]])
            WISFSS.append(key + [WISFS[-1, -1]])
            # scores at every horizon (ranked model rankx)
            scoresByHorizon.append(
                np.column_stack([np.full(forecastingperiod, rankx), MAEFS_model1[:, 0], MAEFS_model1[:, 1],
                                 MSEFS_model1[:, 1], PIFS_model1[:, 1], WISFS[:, 1]])
            )  # fmt: skip

        # <=============================================================================================>
        # <============= Construct ensemble model with the top-ranked models 1 to rankx ================>
        # <=============================================================================================>
        if rankx > 1:
            (RMSECS_model1, MSECS_model1, MAECS_model1, PICS_model1, MISCS_model1, WISC, RMSEFS_model1,
             MSEFS_model1, MAEFS_model1, PIFS_model1, MISFS_model1, WISFS, forecast1, quantilesc,
             quantilesf) = getensemblesubepidemicmodels(
                cumulative1, cadfilename2, datevecfirst1, npatches_fixed, onset_fixed, smoothfactor1, outbreakx,
                cadregion, caddate1, caddisease, datatype, flag1, method1, dist1, calibrationperiod1,
                list(range(1, rankx + 1)), forecastingperiod, getperformance, weight_type1, WISC_hash, WISC_hash,
                printscreen1, ensemblePanelGrid, input_dir=input_dir, output_dir=output_dir,
            )  # fmt: skip

            key = [100 + rankx, outbreakx, dnum]
            RMSECSS.append(key + [RMSECS_model1])
            MSECSS.append(key + [MSECS_model1])
            MAECSS.append(key + [MAECS_model1])
            PICSS.append(key + [PICS_model1])
            MISCSS.append(key + [MISCS_model1])
            WISCSS.append(key + [WISC])

            if scored:
                RMSEFSS.append(key + [RMSEFS_model1[-1, -1]])
                MSEFSS.append(key + [MSEFS_model1[-1, -1]])
                MAEFSS.append(key + [MAEFS_model1[-1, -1]])
                PIFSS.append(key + [PIFS_model1[-1, -1]])
                MISFSS.append(key + [MISFS_model1[-1, -1]])
                WISFSS.append(key + [WISFS[-1, -1]])
                # scores at every horizon (Ensemble(rankx))
                scoresByHorizon.append(
                    np.column_stack([np.full(forecastingperiod, 100 + rankx), MAEFS_model1[:, 0],
                                     MAEFS_model1[:, 1], MSEFS_model1[:, 1], PIFS_model1[:, 1], WISFS[:, 1]])
                )  # fmt: skip

    if deletetempfiles:  # delete the Forecast..mat files after use
        for j in topmodels1:
            try:
                os.remove(names.forecast_mat(forecastingperiod, j))
            except OSError:
                pass

    if printscreen1:
        for fig in (fig10, fig11, plt.figure(500) if len(topmodels1) > 1 else None):
            if fig is not None:
                try:
                    fig.tight_layout()
                except Exception:
                    pass

    def table(rows_):
        return np.array(rows_, dtype=float).reshape(-1, 4)

    MAECSS, MSECSS, PICSS, WISCSS = table(MAECSS), table(MSECSS), table(PICSS), table(WISCSS)
    MAEFSS, MSEFSS, PIFSS, WISFSS = table(MAEFSS), table(MSEFSS), table(PIFSS), table(WISFSS)
    AICc_rank1 = np.array(AICc_rank1, dtype=float).reshape(-1, 2)
    relativelik_rank1 = np.array(relativelik_rank1, dtype=float).reshape(-1, 2)
    ranks = np.array(topmodels1, dtype=float)

    # <=============================================================================================>
    # <=================== Plot forecasting performance metrics of the top-ranked models ===========>
    # <=============================================================================================>
    if scored:
        index1 = np.flatnonzero(MAEFSS[:, 0] >= 100)
        index2 = np.flatnonzero(MAEFSS[:, 0] < 100)

        if printscreen1:
            _score_figure(400, MAEFSS[index2, 0], [MAEFSS[index2, 3], MSEFSS[index2, 3], PIFSS[index2, 3],
                                                    WISFSS[index2, 3]], r"$i_{th}$Ranked Model")  # fmt: skip

        # <=========================================================================================>
        # <============= Save file with top-ranked models' forecasting performance metrics =========>
        # <=========================================================================================>
        performanceTop = np.column_stack([ranks, MAEFSS[index2, 3], MSEFSS[index2, 3], PIFSS[index2, 3],
                                          WISFSS[index2, 3], AICc_rank1[:, 1], relativelik_rank1[:, 1]])  # fmt: skip
        writetable(names.csv("performance-forecasting-topRanked", forecastingperiod), performanceTop,
                   ["i_th-ranked model", "MAE", "MSE", "Coverage 95%PI", "WIS", "AICc", "RelativeLikelihood"])  # fmt: skip

        sizes = 1 + np.arange(1, index1.size + 1, dtype=float)
        if index1.size and printscreen1:
            _score_figure(401, sizes, [MAEFSS[index1, 3], MSEFSS[index1, 3], PIFSS[index1, 3], WISFSS[index1, 3]],
                          "Ensemble(i) model")  # fmt: skip

        # <=========================================================================================>
        # <================= Save file with ensemble forecasting performance metrics ===============>
        # <=========================================================================================>
        performanceEns = np.column_stack([sizes, MAEFSS[index1, 3], MSEFSS[index1, 3], PIFSS[index1, 3],
                                          WISFSS[index1, 3]])  # fmt: skip
        writetable(names.csv("performance-forecasting-Ensemble", forecastingperiod, weight_type1), performanceEns,
                   ["Ensemble(i) model", "MAE", "MSE", "Coverage 95%PI", "WIS"])  # fmt: skip

        # <=========================================================================================>
        # <=================== Baseline forecast, and scores at every horizon ======================>
        # <=========================================================================================>
        # A model forecast means little without a reference. The baseline is the
        # mean of the last observations, repeated over the forecast period, with
        # an interval from the errors that rule made earlier in the same series.
        baselineWindow = fitSettings.baselineWindow
        if baselineWindow is None:
            baselineWindow = 7 if DT == 1 else 1  # daily data: one week; otherwise the last observation
        baseline = subepidemicBaselineForecast(data1[:, 1], forecastingperiod, baselineWindow)
        baselineCurves1 = np.concatenate([data1[:, 1], baseline.point]).reshape(-1, 1)
        baselineCurves2 = np.vstack([np.tile(data1[:, 1].reshape(-1, 1), (1, baseline.samples.shape[1])),
                                     baseline.samples])  # fmt: skip
        r = computeforecastperformance(data1, datalatest2, baselineCurves1, baselineCurves2, forecastingperiod)
        MSEFS_base, MAEFS_base, PIFS_base = r[6], r[7], r[8].copy()
        if np.all(np.isfinite(baseline.samples)):
            _, WISFS_base, _, _ = computeWIS(data1, datalatest2, baselineCurves2, forecastingperiod)
        else:
            WISFS_base = np.column_stack([np.arange(1, forecastingperiod + 1), np.full(forecastingperiod, np.nan)])
            PIFS_base[:, 1] = np.nan
            subepidemicWarning(
                "SubEpiPredict:BaselineWithoutInterval",
                "The calibration series is too short to give the flat baseline an interval at every horizon "
                "(%d observations, window %d, horizon %d). Its MAE and MSE are reported; its WIS and coverage "
                "are not.",
                data1.shape[0], baseline.window, forecastingperiod,
            )  # fmt: skip
        baseRows = np.column_stack([np.zeros(forecastingperiod), MAEFS_base[:, 0], MAEFS_base[:, 1],
                                    MSEFS_base[:, 1], PIFS_base[:, 1], WISFS_base[:, 1]])  # fmt: skip
        scoresAll = np.vstack([baseRows] + scoresByHorizon)
        scoresAll = scoresAll[np.lexsort((scoresAll[:, 1], scoresAll[:, 0])), :]  # sortrows(...,[1 2])
        with np.errstate(divide="ignore", invalid="ignore"):
            relativeWIS = scoresAll[:, 5] / WISFS_base[scoresAll[:, 1].astype(int) - 1, 1]
        scoresByHorizonOut = np.column_stack([scoresAll, relativeWIS])
        writetable(names.csv("performance-forecasting-byHorizon", forecastingperiod, weight_type1),
                   scoresByHorizonOut,
                   ["model", "horizon", "MAE", "MSE", "Coverage 95%PI", "WIS", "WIS relative to baseline"])  # fmt: skip
        performanceBaseline = np.array([MAEFS_base[-1, 1], MSEFS_base[-1, 1], PIFS_base[-1, 1], WISFS_base[-1, 1]])

        if verbose:
            rowsEnd = scoresByHorizonOut[scoresByHorizonOut[:, 1] == forecastingperiod, :]
            print(
                "\nForecast scores over %d time unit(s). Baseline: flat at %.6g, the mean of the last %d observation(s)."
                % (forecastingperiod, baseline.level, baseline.window)
            )
            print("  %-14s %10s %10s %10s %14s" % ("model", "MAE", "coverage", "WIS", "WIS/baseline"))
            for row in rowsEnd:
                if row[0] == 0:
                    modelLabel = "baseline"
                elif row[0] > 100:
                    modelLabel = "Ensemble(%d)" % (row[0] - 100)
                else:
                    modelLabel = "ranked %d" % row[0]
                print("  %-14s %10.2f %10.1f %10.2f %14.2f" % (modelLabel, row[2], row[4], row[5], row[6]))
            print(
                "  (scores of every horizon from 1 to %d are in the file performance-forecasting-byHorizon-... in %s)"
                % (forecastingperiod, output_dir)
            )
    elif getperformance == 1:
        performanceTop = np.zeros((0, 7))

    # <=============================================================================================>
    # <============= Save file with top-ranked models' calibration performance metrics =============>
    # <=============================================================================================>
    index2 = np.flatnonzero(MAECSS[:, 0] < 100)
    performanceC = np.column_stack([ranks, MAECSS[index2, 3], MSECSS[index2, 3], PICSS[index2, 3], WISCSS[index2, 3]])
    writetable(names.csv("performance-calibration-topRanked"), performanceC,
               ["i_th-ranked model", "MAE", "MSE", "Coverage 95%PI", "WIS"])  # fmt: skip

    # <=============================================================================================>
    # <================= Save file with ensemble calibration performance metrics ===================>
    # <=============================================================================================>
    index1 = np.flatnonzero(MAECSS[:, 0] > 100)
    performanceEnsC = np.column_stack([MAECSS[index1, 0] - 100, MAECSS[index1, 3], MSECSS[index1, 3],
                                       PICSS[index1, 3], WISCSS[index1, 3]])  # fmt: skip
    writetable(names.csv("performance-calibration-Ensemble", forecastingperiod, weight_type1), performanceEnsC,
               ["Ensemble(i) model", "MAE", "MSE", "Coverage 95%PI", "WIS"])  # fmt: skip

    # List the toolbox's warnings of this run once more.
    if verbose:
        subepidemicWarningLog("summary", "plotForecast_subepidemicFramework")

    return performanceTop, performanceEns, performanceBaseline, scoresByHorizonOut


def _score_figure(number, x, columns, xlabel):
    """2-by-2 figure of MAE, MSE, coverage and WIS against the model."""
    import matplotlib.pyplot as plt

    fig = plt.figure(number, figsize=(10, 7.5))
    fig.clf()
    for k, (values, label) in enumerate(zip(columns, ["MAE", "MSE", "Coverage of the 95% PI", "WIS"])):
        ax = fig.add_subplot(2, 2, k + 1)
        ax.plot(x, values, "k-o", linewidth=2)
        ax.set_xlabel(xlabel, fontsize=13)
        ax.set_ylabel(label, fontsize=13)
        ax.tick_params(labelsize=11)
    fig.set_facecolor("white")
    fig.tight_layout()
