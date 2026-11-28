"""Plots and tables of a saved fit: model fit, rankings, reproduction number,
and the standalone curve plotters.

Ports of ``plotFit_subepidemicFramework.m``,
``plotRankings_subepidemicFramework.m``, ``plotReproductionNumber.m``,
``plot_nsubepidemic.m`` and ``plot_nsubepidemic_profiles.m``.

Figures are numbered as in MATLAB and are never shown by the functions: call
``matplotlib.pyplot.show()`` or save them (``plt.figure(301).savefig(...)``).
"""

from __future__ import annotations

import numpy as np

from .errors import SubEpiPredictError, subepidemicWarning, subepidemicWarningLog
from .figures import (
    COMPONENT_COLORS,
    GRAY7,
    RANKING_COLORS,
    GetAdjustedFontSize,
    date_ticks,
    histogram,
    mathLabel,
    plot_curves,
    subepidemicEstimateLine,
    subepidemicPanelGrid,
    subepidemicParameterLabel,
)
from .forecast import (
    DOUBLING_NAMES,
    _doubling_summary,
    _interval_panel,
    load_ranked,
    subepidemicAvailableTopModels,
    subepidemicSkippedRefits,
)
from .identifiability import sci_interval, subepidemicIdentifiability, subepidemicPrintPII, subepidemicWritePII
from .io import OutputNames, load_mat, save_mat, scalar, vector, writetable
from .model import simulateSubepidemic, subepidemicSimulateRefit
from .noise import AddPoissonError
from .optionsets import options as default_options
from .optionsets import options_forecast as default_options_forecast
from .optionsets import options_Rt as default_options_Rt
from .optionsets import resolve_options
from .rt import get_Rt
from .scores import _clip0, computeforecastperformance, computeWIS
from .search import subepidemicSearchSpace
from .search import unique_rows_stable as unique_rows
from .settings import resolve_fit_settings, subepidemicSeed
from .utils import num2ordinal, num2str, quantile


def _area_and_date(opts, outbreakx_pass, caddate1_pass):
    outbreakx = opts.outbreakx if outbreakx_pass is None or np.size(outbreakx_pass) == 0 else outbreakx_pass
    outbreakx = int(np.asarray(outbreakx).ravel()[0])
    caddate1 = opts.caddate1 if caddate1_pass is None or caddate1_pass == "" else caddate1_pass
    return outbreakx, caddate1


def _ci(values):
    """[mean, 2.5% quantile, 97.5% quantile] of a sample."""
    return [float(np.mean(values)), quantile(values, 0.025), quantile(values, 0.975)]


# <============================================================================>
# <================================ Model fit =================================>
# <============================================================================>


def plotFit_subepidemicFramework(
    outbreakx_pass=None,
    caddate1_pass=None,
    options=None,
    fit_settings=None,
    output_dir="./output",
    plot=True,
    verbose=True,
):
    """Plot the fit of the top-ranked models and write their parameter tables.

    Port of ``plotFit_subepidemicFramework.m``. Run
    :func:`Run_Fit_subepidemicFramework` with the same options first.

    Figures: 300+k, the bootstrap distributions of r, p, a and K of every
    sub-epidemic of ranked model k (red dashed line: 95% interval; black
    dashed line: the estimate); 300+n+1 for n ranked models, one row per
    model with the fit and bootstrap data, the sub-epidemic components and
    the residuals; 400, the calibration scores.

    Files (names as in MATLAB): ``param-r-``, ``param-p-``, ``param-a-``,
    ``param-K-``, ``param-NB-alpha-``, ``param-NB-d-``, ``PII-`` and
    ``doublingTimes-ranked(k)-...csv`` for each ranked model,
    ``performance-calibration-topRanked-`` and ``PII-topRanked-...csv``.

    Returns ``performance``: one row per ranked model, ``[rank, MAE, MSE,
    coverage of the 95% PI, WIS, AICc, relative likelihood]`` over the
    calibration period.
    """
    opts = resolve_options(options, default_options)
    fitSettings = resolve_fit_settings(fit_settings)
    outbreakx, caddate1 = _area_and_date(opts, outbreakx_pass, caddate1_pass)

    caddisease = opts.caddisease
    datatype = opts.datatype
    DT = opts.DT
    method1 = opts.method1
    onset_fixed = opts.onset_fixed

    names = OutputNames(opts, outbreakx, caddate1, output_dir)

    topmodels1 = list(range(1, opts.topmodelsx + 1))
    if opts.npatches_fixed == 1:
        topmodels1 = [1]

    # A fresh list of warnings for this run.
    subepidemicWarningLog("start")
    topmodels1 = subepidemicAvailableTopModels(topmodels1, names)

    # Random seed (setting rngSeed), for anything below that draws random numbers.
    subepidemicSeed("plotFit_subepidemicFramework", fitSettings, verbose=verbose)

    MAES, MSES, PIS, WISS = [], [], [], []
    AICc_rank1 = []
    relativelik_rank1 = []
    piiRanked = []  # PII metrics of each ranked model

    if plot:
        import matplotlib.pyplot as plt

        plt.close("all")
        fs = GetAdjustedFontSize()
        nModels = len(topmodels1)
        figFit = plt.figure(300 + topmodels1[-1] + 1, figsize=(19, 4.3 * nModels))

    for cc, rank1 in enumerate(topmodels1):

        S = load_ranked(names, rank1)
        npatches = S.npatches
        onset_thr = S.onset_thr
        flag1 = S.flag1
        Phatss = S.Phatss
        data1 = S.data1
        data = S.data
        curves = S.curves
        M = S.M

        # Residuals are the unchanged saved calibration fit minus the observed data.
        bestfit = S.bestfit
        if bestfit is None or bestfit.size != data1.shape[0] or np.any(~np.isfinite(bestfit)):
            raise SubEpiPredictError(
                "plotFit_subepidemicFramework:InvalidSavedBestFit",
                "Rank %d must contain a finite real bestfit vector with one value per calibration data row." % rank1,
            )
        bestfit_original = bestfit

        # The fitted estimate is the ranked vector P0. It is the number shown in
        # every title; the mean of the bootstrap refits is not used for it.
        if S.P0 is None or S.P0.size != 4 * npatches + 2 or np.any(~np.isfinite(S.P0)):
            raise SubEpiPredictError(
                "plotFit_subepidemicFramework:InvalidSavedEstimate",
                "Rank %d must contain the fitted parameter vector P0 with %d finite values. Rerun the fit."
                % (rank1, 4 * npatches + 2),
            )
        fittedEstimate = S.P0

        AICc_rank1.append([rank1, S.AICc_best])
        relativelik_rank1.append([rank1, S.relativelik_i[rank1 - 1]])

        timevect = data1[:, 0]

        rs = Phatss[:, 0:npatches]
        ps = Phatss[:, npatches : 2 * npatches]
        as_ = Phatss[:, 2 * npatches : 3 * npatches]
        Ks = Phatss[:, 3 * npatches : 4 * npatches]

        param_alpha = np.array([_ci(Phatss[:, -2])])
        param_d = np.array([_ci(Phatss[:, -1])])

        param_r = np.array([_ci(rs[:, j]) for j in range(npatches)])
        param_p = np.array([_ci(ps[:, j]) for j in range(npatches)])
        param_a = np.array([_ci(as_[:, j]) for j in range(npatches)])
        param_K = np.array([_ci(Ks[:, j]) for j in range(npatches)])

        if verbose:
            print(subepidemicParameterLabel("alpha", fittedEstimate[-2], param_alpha[0, 1], param_alpha[0, 2]))
            print(subepidemicParameterLabel("d", fittedEstimate[-1], param_d[0, 1], param_d[0, 2]))

        # <========================================================================================>
        # <======================= Plot parameter estimates =======================================>
        # <========================================================================================>
        if plot:
            figPar = plt.figure(300 + rank1, figsize=(19, 3.9 * npatches))
        for j in range(npatches):
            blocks = (
                ("r", rs[:, j], param_r[j], fittedEstimate[j]),
                ("p", ps[:, j], param_p[j], fittedEstimate[npatches + j]),
                ("a", as_[:, j], param_a[j], fittedEstimate[2 * npatches + j]),
                ("K", Ks[:, j], param_K[j], fittedEstimate[3 * npatches + j]),
            )
            for k, (letter, values, ci, estimate) in enumerate(blocks):
                label = subepidemicParameterLabel("%s_%d" % (letter, j + 1), estimate, ci[1], ci[2])
                if verbose and letter != "r":
                    print(label)
                if not plot:
                    continue
                ax = figPar.add_subplot(npatches, 4, 4 * j + k + 1)
                histogram(ax, values)
                ax.set_xlabel(letter, fontsize=fs)
                ax.set_title(mathLabel(label), fontsize=fs, fontweight="bold")
                ax.plot([ci[1], ci[2]], [10, 10], "r--", linewidth=2)
                subepidemicEstimateLine(float(estimate), ax)
                ax.tick_params(labelsize=fs - 2)
                if k == 0:
                    from matplotlib.patches import Patch

                    ax.legend(handles=[Patch(facecolor=(0.24, 0.15, 0.66), edgecolor="k")],
                              labels=["%s Ranked Model" % num2ordinal(rank1)], fontsize=fs - 3)  # fmt: skip
        if plot:
            figPar.set_facecolor("white")
            figPar.tight_layout()

        # <========================================================================================>
        # <=========================== Save file with parameters ==================================>
        # <========================================================================================>
        # Columns: mean and 95% limits of the bootstrap refits, the width of the
        # interval on the log10 scale (SCI), the estimate (the value of the ranked
        # fit, the centre of the bootstrap) and the median of the refits.
        forecastingperiod = 0
        sub = np.arange(1, npatches + 1, dtype=float)
        for letter, block, col0 in (("r", param_r, 0), ("p", param_p, npatches), ("a", param_a, 2 * npatches),
                                    ("K", param_K, 3 * npatches)):  # fmt: skip
            tableValues = np.column_stack([
                sub, block, sci_interval(block[:, 1], block[:, 2]), fittedEstimate[col0 : col0 + npatches],
                np.median(Phatss[:, col0 : col0 + npatches], axis=0),
            ])  # fmt: skip
            writetable(
                names.csv("param-%s-ranked(%s)" % (letter, num2str(rank1)), forecastingperiod), tableValues,
                ["sub-epidemic", "%s mean" % letter, "%s 95%% CI LB" % letter, "%s 95%% CI UB" % letter,
                 "%s SCI" % letter, "%s estimate" % letter, "%s median" % letter],
            )  # fmt: skip

        def errorTable(letter, block, column):
            writetable(
                names.csv("param-NB-%s-ranked(%s)" % (letter, num2str(rank1)), forecastingperiod),
                np.column_stack([block, sci_interval(block[:, 1], block[:, 2]), fittedEstimate[column],
                                 np.median(Phatss[:, column])]),
                ["%s mean" % letter, "%s 95%% CI LB" % letter, "%s 95%% CI UB" % letter, "%s SCI" % letter,
                 "%s estimate" % letter, "%s median" % letter],
            )  # fmt: skip

        if method1 in (3, 4, 5):  # save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2
            errorTable("alpha", param_alpha, -2)
        if method1 == 5:
            errorTable("d", param_d, -1)

        # <========================================================================================>
        # <======================= Parameter identifiability (PII) ================================>
        # <========================================================================================>
        # The search bounds are those saved with the fit; a fit saved without
        # them gets the bounds of the current settings, with a warning.
        diag = S.get("bootstrapDiagnostics")
        nz = 4 * npatches + 2
        if isinstance(diag, dict) and np.size(diag.get("searchLB", [])) == nz and np.size(diag.get("searchUB", [])) == nz:
            piiLB = vector(diag["searchLB"])
            piiUB = vector(diag["searchUB"])
        else:
            piiLB, piiUB, _, _ = subepidemicSearchSpace(npatches, data1[:, 1], flag1, method1, fitSettings)
            subepidemicWarning(
                "SubEpiPredict:PIIBoundsRebuilt",
                "Ranked model %d was saved without its search bounds. The PII report uses the bounds of the "
                "current settings; if the bound settings changed since the fit, rerun the fit.",
                rank1,
            )
        pii = subepidemicIdentifiability(Phatss, piiLB, piiUB, npatches, fittedEstimate, fitSettings.PIIscale)
        piiRanked.append(pii)
        if verbose:
            subepidemicPrintPII(pii, "Parameter identifiability (PII) of ranked model %d" % rank1)
        subepidemicWritePII(names.csv("PII-ranked(%s)" % num2str(rank1), forecastingperiod), pii)

        # <========================================================================================>
        # <================================ Plot model fit ========================================>
        # <========================================================================================>
        LB1 = quantile(curves, 0.025, axis=1)
        LB1 = _clip0(LB1)
        UB1 = quantile(curves, 0.975, axis=1)
        UB1 = _clip0(UB1)
        top = float(np.max(data[:, 1])) * 1.3

        fittedcurves = np.zeros((timevect.size, M))
        components = [[] for _ in range(npatches)]
        for realization in range(M):
            _, x, totinc, patchIncidence, _ = simulateSubepidemic(
                timevect, data1[0, 1], npatches, onset_fixed, onset_thr, flag1,
                Phatss[realization, 0:npatches], Phatss[realization, npatches : 2 * npatches],
                Phatss[realization, 2 * npatches : 3 * npatches], Phatss[realization, 3 * npatches : 4 * npatches],
            )  # fmt: skip
            for j in range(npatches):
                components[j].append(patchIncidence[:, j])
            fittedcurves[:, realization] = totinc

        resid1 = bestfit_original - data1[:, 1]

        if plot:
            ax = figFit.add_subplot(nModels, 3, 3 * cc + 1)
            plot_curves(ax, timevect, curves, "c")
            ax.plot(timevect, np.median(curves, axis=1), "r-", linewidth=2)
            ax.plot(timevect, LB1, "r--", linewidth=2)
            ax.plot(timevect, UB1, "r--", linewidth=2)
            ax.plot(timevect, data1[:, 1], "ko", markerfacecolor="none", markeredgewidth=1.5)
            ax.set_xlabel("Time", fontsize=fs)
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.set_xlim(timevect[0], timevect[-1] + 1)
            ax.set_ylim(0, top)
            ax.tick_params(labelsize=fs - 2)

            ax = figFit.add_subplot(nModels, 3, 3 * cc + 2)
            for j in range(npatches):
                lc = plot_curves(ax, timevect, np.column_stack(components[j]), COMPONENT_COLORS[j % 7])
                if j == 0 and lc is not None:
                    lc.set_label("Sub-epidemics=%s; $C_{thr}$=%s" % (num2str(npatches), num2str(onset_thr)))
            plot_curves(ax, timevect, fittedcurves, GRAY7)
            ax.set_xlabel("Time (days)", fontsize=fs)
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.plot(data[:, 0], data[:, 1], "ko", markerfacecolor="none", markeredgewidth=1.5)
            ax.set_xlim(timevect[0], timevect[-1] + 1)
            ax.set_ylim(0, top)
            ax.legend(fontsize=fs - 3, loc="upper left")
            ax.set_title("%s Ranked Model" % num2ordinal(rank1), fontsize=fs, fontweight="bold")
            ax.tick_params(labelsize=fs - 2)

            ax = figFit.add_subplot(nModels, 3, 3 * cc + 3)
            markerline, stemlines, baseline_ = ax.stem(timevect, resid1, linefmt="b-", markerfmt="bo", basefmt="k-")
            markerline.set_markerfacecolor("none")
            ax.set_xlabel("Time (days)", fontsize=fs)
            ax.set_ylabel("Residuals (fitted - observed)", fontsize=fs)
            ax.set_xlim(timevect[0], timevect[-1] + 1)
            ax.set_ylim(float(np.min(resid1)) - 1, float(np.max(resid1)) + 1)
            ax.tick_params(labelsize=fs - 2)

        # compute doubling times
        seq_doublingtimes, param_doubling = _doubling_summary(fittedcurves, DT, M)
        if verbose:
            print("Doubling time of ranked model %d, mean (95%% CI): %.4g (%.4g, %.4g)"
                  % (rank1, param_doubling[0], param_doubling[1], param_doubling[2]))  # fmt: skip

        # <=============================================================================================>
        # <============================== Save file with doubling time estimates =======================>
        # <=============================================================================================>
        writetable(names.csv("doublingTimes-ranked(%s)" % num2str(rank1), forecastingperiod), seq_doublingtimes,
                   DOUBLING_NAMES)  # fmt: skip

        r = computeforecastperformance(data1, data1, fittedcurves, curves, 0)
        WISC, _, _, _ = computeWIS(data1, data1, curves, 0)

        MAES.append(r[2])
        MSES.append(r[1])
        PIS.append(r[3])
        WISS.append(WISC)

    AICc_rank1 = np.array(AICc_rank1, dtype=float).reshape(-1, 2)
    relativelik_rank1 = np.array(relativelik_rank1, dtype=float).reshape(-1, 2)
    ranks = np.array(topmodels1, dtype=float)

    if plot:
        figFit.set_facecolor("white")
        figFit.tight_layout()
        from .forecast import _score_figure

        _score_figure(400, ranks, [MAES, MSES, PIS, WISS], r"$i_{th}$Ranked Model")

    # <=============================================================================================>
    # <============= Save file with top-ranked models' performance metrics (calibration) ===========>
    # <=============================================================================================>
    performance = np.column_stack([ranks, MAES, MSES, PIS, WISS, AICc_rank1[:, 1], relativelik_rank1[:, 1]])
    writetable(names.csv("performance-calibration-topRanked"), performance,
               ["i_th-ranked model", "MAE", "MSE", "Coverage 95%PI", "WIS", "AICc", "RelativeLikelihood"])  # fmt: skip

    # One file with the PII metrics of every ranked model, for comparing models.
    subepidemicWritePII(names.csv("PII-topRanked"), piiRanked, topmodels1)

    # List the toolbox's warnings of this run once more.
    if verbose:
        subepidemicWarningLog("summary", "plotFit_subepidemicFramework")

    return performance


# <============================================================================>
# <================================= Rankings =================================>
# <============================================================================>


def plotRankings_subepidemicFramework(
    outbreakx_pass=None,
    caddate1_pass=None,
    options=None,
    output_dir="./output",
    plot=True,
    verbose=True,
):
    """Plot the AICc ranking and the sub-epidemic profiles of the top-ranked models.

    Port of ``plotRankings_subepidemicFramework.m``. Reads the ranking
    (``ABC-ensem-...mat``) saved by the fit.

    Figures: 99, AICc, relative likelihood and evidence ratio of the ranked
    models; 100, one panel per ranked model with its sub-epidemic components
    (colors), their sum (black) and the data.

    Files: ``AICc-topRanked-...csv`` and ``AICc_bests-...mat``.

    Returns ``performanceTop``: one row per ranked model, ``[rank, AICc,
    relative likelihood]``. (The MATLAB function returns nothing.)
    """
    opts = resolve_options(options, default_options)
    outbreakx, caddate1 = _area_and_date(opts, outbreakx_pass, caddate1_pass)

    # A fresh list of warnings for this run.
    subepidemicWarningLog("start")

    caddisease = opts.caddisease
    datatype = opts.datatype
    onset_fixed = opts.onset_fixed
    flag1 = opts.flag1
    names = OutputNames(opts, outbreakx, caddate1, output_dir)

    topmodelsx = list(range(1, opts.topmodelsx + 1))

    abc = load_mat(names.abc_mat())
    RMSES = np.asarray(abc.RMSES, dtype=float).reshape(-1, 3)
    PS = np.asarray(abc.PS, dtype=float)
    relativelik_i = vector(abc.relativelik_i)
    timevect = vector(abc.timevect)
    I0 = scalar(abc.I0)
    data = vector(abc.data)

    # remove repeated rows
    RMSES, index1 = unique_rows(RMSES)
    PS = PS[index1, :]

    # The options can ask for more top models than the screen retained.
    if topmodelsx[-1] > RMSES.shape[0]:
        subepidemicWarning(
            "SubEpiPredict:TopModelsCapped",
            "The options ask for %d top-ranked models, but the screen retained %d. Ranks 1 to %d are used.",
            topmodelsx[-1], RMSES.shape[0], RMSES.shape[0],
        )  # fmt: skip
        topmodelsx = list(range(1, RMSES.shape[0] + 1))
    idx = np.array(topmodelsx) - 1
    ranks = np.array(topmodelsx, dtype=float)

    AICc_bests = RMSES[idx, 2]

    if plot:
        import matplotlib.pyplot as plt

        plt.close("all")
        fs = GetAdjustedFontSize()

        # plot relative likelihoods of the models
        fig = plt.figure(99, figsize=(15, 4.5))
        for k, (values, label) in enumerate(
            ((RMSES[idx, 2], "AICc"), (relativelik_i[idx], "Relative likelihood"),
             (1.0 / relativelik_i[idx], "Evidence ratio"))
        ):  # fmt: skip
            ax = fig.add_subplot(1, 3, k + 1)
            ax.plot(ranks, values, "ko-", linewidth=2)
            ax.set_xlabel(r"$i_{th}$ Ranked Model", fontsize=fs)
            ax.set_ylabel(label, fontsize=fs)
            ax.tick_params(labelsize=fs - 2)
        fig.set_facecolor("white")
        fig.tight_layout()

        rows, cols = subepidemicPanelGrid(len(topmodelsx))
        fig100 = plt.figure(100, figsize=(7.5 * cols, 4.2 * rows))

    for rank in topmodelsx:
        npatches = int(RMSES[rank - 1, 0])
        onset_thr = float(RMSES[rank - 1, 1])
        AICc_best = float(RMSES[rank - 1, 2])
        P = PS[rank - 1, : npatches * 4 + 2]

        rs_hat = P[0:npatches]
        ps_hat = P[npatches : 2 * npatches]
        as_hat = P[2 * npatches : 3 * npatches]
        Ks_hat = P[3 * npatches : 4 * npatches]

        _, x, totinc, patchIncidence, _ = simulateSubepidemic(
            timevect, I0, npatches, onset_fixed, onset_thr, flag1, rs_hat, ps_hat, as_hat, Ks_hat
        )

        if verbose:
            print("Ranked model %d: %d sub-epidemic(s), C_thr = %s, AICc = %s"
                  % (rank, npatches, num2str(onset_thr), num2str(AICc_best, 10)))  # fmt: skip
            print("   r = %s\n   p = %s\n   a = %s\n   K = %s"
                  % tuple(" ".join("%.5g" % v for v in block) for block in (rs_hat, ps_hat, as_hat, Ks_hat)))  # fmt: skip

        if plot:
            ax = fig100.add_subplot(rows, cols, rank)
            for j in range(npatches):
                (line1,) = ax.plot(timevect, patchIncidence[:, j], RANKING_COLORS[j % 6] + "-", linewidth=2)
                if j == 0:
                    line1.set_label("Sub-epidemics=%s; $C_{thr}$=%s" % (num2str(npatches), num2str(onset_thr)))
            ax.plot(timevect, data, "bo", markerfacecolor="none", markeredgewidth=1.5)
            ax.plot(timevect, totinc, "k", linewidth=3)
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.set_title("%s  Ranked Model; AICc=%s" % (num2ordinal(rank), num2str(AICc_best, 6)), fontsize=fs,
                         fontweight="bold")  # fmt: skip
            ax.legend(fontsize=fs - 3)
            ax.tick_params(labelsize=fs - 2)
            if rank > len(topmodelsx) - cols:
                ax.set_xlabel("Time", fontsize=fs)

    if plot:
        fig100.set_facecolor("white")
        fig100.tight_layout()

    save_mat(names.aicc_bests_mat(topmodelsx[-1]), dict(AICc_bests=AICc_bests.reshape(-1, 1)))

    # <=============================================================================================>
    # <============= Save file with top-ranked models' AICc and relative likelihood ================>
    # <=============================================================================================>
    performanceTop = np.column_stack([ranks, RMSES[idx, 2], relativelik_i[idx]])
    writetable(names.csv("AICc-topRanked"), performanceTop, ["i_th-ranked model", "AICc", "RelativeLikelihood"])

    # List the toolbox's warnings of this run once more.
    if verbose:
        subepidemicWarningLog("summary", "plotRankings_subepidemicFramework")

    return performanceTop


# <============================================================================>
# <============================ Reproduction number ===========================>
# <============================================================================>


def plotReproductionNumber(
    outbreakx_pass=None,
    caddate1_pass=None,
    options=None,
    options_forecast=None,
    options_Rt=None,
    fit_settings=None,
    output_dir="./output",
    plot=True,
    verbose=True,
):
    """Effective reproduction number of the top-ranked models.

    Port of the script ``plotReproductionNumber.m``. For every bootstrap
    parameter vector of a ranked model the curve is simulated over the
    calibration period and the forecast horizon of ``options_forecast``, and
    Rt is computed from it with the renewal equation and the generation
    interval of ``options_Rt`` (:func:`get_Rt`).

    Figures: 10, trajectories and components; 502, Rt (cyan: every curve;
    red: median and 95% limits); 11, the forecast with its 95% interval.

    Files: ``Rt-ranked(k)-...csv`` with the columns time, Rt median, 95% CI
    LB and UB at the whole time units.

    Returns ``param_Rt``: one row per ranked model, ``[last time, median, 95%
    LB, 95% UB, standard deviation]`` of Rt at the end of the curve. (The
    MATLAB script returns nothing; the area and the date are those of the
    options there.)
    """
    opts = resolve_options(options, default_options)
    fopts = resolve_options(options_forecast, default_options_forecast)
    ropts = resolve_options(options_Rt, default_options_Rt)
    fitSettings = resolve_fit_settings(fit_settings)
    outbreakx, caddate1 = _area_and_date(opts, outbreakx_pass, caddate1_pass)

    caddisease = opts.caddisease
    datatype = opts.datatype
    datevecfirst1 = opts.datevecfirst1
    DT = opts.DT
    method1 = opts.method1
    dist1 = opts.dist1
    onset_fixed = opts.onset_fixed

    names = OutputNames(opts, outbreakx, caddate1, output_dir)

    topmodels1 = list(range(1, opts.topmodelsx + 1))
    if opts.npatches_fixed == 1:
        topmodels1 = [1]

    # A fresh list of warnings for this run.
    subepidemicWarningLog("start")
    topmodels1 = subepidemicAvailableTopModels(topmodels1, names)

    # Random seed (setting rngSeed): the forecast noise repeats from run to run.
    subepidemicSeed("plotReproductionNumber", fitSettings, verbose=verbose)

    forecastingperiod = int(fopts.forecastingperiod)
    type_GId1 = ropts.type_GId1
    mean_GI1 = ropts.mean_GI1
    var_GI1 = ropts.var_GI1
    minGIcoverage = ropts.minGIcoverage

    param_Rt = []
    nModels = len(topmodels1)

    if plot:
        import matplotlib.pyplot as plt

        plt.close("all")
        fs = GetAdjustedFontSize()
        fig10 = plt.figure(10, figsize=(6.2 * nModels, 4.6))
        fig502 = plt.figure(502, figsize=(6.2 * nModels, 4.6))
        fig11 = plt.figure(11, figsize=(6.2 * nModels, 4.6))

    for cc1, rankx in enumerate(topmodels1, start=1):

        S = load_ranked(names, rankx)
        npatches = S.npatches
        onset_thr = S.onset_thr
        flag1 = S.flag1
        Phatss = S.Phatss
        data1 = S.data1
        timelags = S.timelags
        rank1 = S.rank1
        factor1 = S.factor1
        M = S.M

        d_hat = 1
        timevect = data1[:, 0]
        timevect2 = np.arange(0, S.t_window[-1] + forecastingperiod, dtype=float)

        cols1, cols2, Rcols = [], [], []
        components = [[] for _ in range(npatches)]
        ts = None

        nSimulationsSkipped = 0  # refits left out because their curve could not be simulated
        firstSkipMessage = ""

        for realization in range(M):
            rs_hat = Phatss[realization, 0:npatches]
            ps_hat = Phatss[realization, npatches : 2 * npatches]
            as_hat = Phatss[realization, 2 * npatches : 3 * npatches]
            Ks_hat = Phatss[realization, 3 * npatches : 4 * npatches]

            # methods 3, 4 and 5 estimate alpha; only method 5 estimates d
            if method1 >= 3 and method1 != 6:
                factor1 = Phatss[realization, -2]
            if method1 == 5:
                d_hat = Phatss[realization, -1]

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
            cols2.append(AddPoissonError(np.cumsum(totinc), 20, dist1, factor1, d_hat))

            Rs, ts, _ = get_Rt(timevect2, totinc, type_GId1, mean_GI1, var_GI1, minGIcoverage)
            Rcols.append(Rs)

        subepidemicSkippedRefits("Rt", rank1, nSimulationsSkipped, M, firstSkipMessage, fitSettings)

        curvesforecasts1 = np.column_stack(cols1)
        curvesforecasts2 = np.hstack(cols2)
        Rss = np.column_stack(Rcols)

        # Rt is blank where the curve covers too little of the generation interval.
        hasValue = np.any(np.isfinite(Rss), axis=1)
        if not np.any(hasValue):
            subepidemicWarning(
                "SubEpiPredict:RtNotEstimable",
                "Ranked model %d: the curve is shorter than the generation interval needs "
                "(coverage below %g everywhere), so no Rt value is reported. Check the units "
                "and values in options_Rt.",
                rank1, minGIcoverage,
            )  # fmt: skip
        elif verbose:
            firstRt = int(np.flatnonzero(hasValue)[0])
            print(
                "Ranked model %d: Rt is reported from time %.2f on (generation-interval coverage of at least %g)."
                % (rank1, ts[firstRt + 1], minGIcoverage)
            )

        Rtmedian = quantile(Rss, 0.5, axis=1)
        RtLB = quantile(Rss, 0.025, axis=1)
        RtUB = quantile(Rss, 0.975, axis=1)

        # <=============================================================================================>
        # <================================ Save file with Rt estimates ================================>
        # <=============================================================================================>
        # One row per whole time unit: Rss[k] is the estimate at time ts[k+1].
        nWhole = int(np.floor(round((ts[-1] - ts[0]) * 100) / 100))
        rowsRt = np.round(np.arange(1, nWhole + 1) * 100).astype(int)
        Rtdata = np.column_stack([ts[rowsRt], Rtmedian[rowsRt - 1], RtLB[rowsRt - 1], RtUB[rowsRt - 1]])
        writetable(names.csv("Rt-ranked(%s)" % num2str(rank1), forecastingperiod), Rtdata,
                   ["time", "Rt median", "Rt 95%CI LB", "Rt 95% CI UB"])  # fmt: skip

        last = Rss[-1, :]
        param_Rt.append([ts[-1], np.median(last), quantile(last, 0.025), quantile(last, 0.975),
                         np.std(last, ddof=1) if last.size > 1 else 0.0])  # fmt: skip

        LB1 = quantile(curvesforecasts2, 0.025, axis=1)
        LB1 = _clip0(LB1)
        UB1 = quantile(curvesforecasts2, 0.975, axis=1)
        UB1 = _clip0(UB1)

        if plot:
            ax = fig10.add_subplot(1, nModels, cc1)
            for j in range(npatches):
                plot_curves(ax, timevect2, np.column_stack(components[j]), COMPONENT_COLORS[j % 7])
            plot_curves(ax, timevect2, curvesforecasts1, GRAY7)
            ax.plot(data1[:, 0], data1[:, 1], "ko", markerfacecolor="none", markeredgewidth=1.5)
            top = float(np.max(data1[:, 1])) * 4
            ax.set_xlim(0, len(timevect2) - 1)
            ax.set_ylim(0, top)
            ax.plot([timevect[-1], timevect[-1]], [0, top], "k--", linewidth=2)
            date_ticks(ax, datevecfirst1, timelags, DT, len(timevect2))
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.set_title("%s Ranked Model" % num2ordinal(rank1), fontsize=fs, fontweight="bold")

            ax = fig502.add_subplot(1, nModels, cc1)
            step = 10  # every tenth grid point is enough for the picture
            plot_curves(ax, ts[1:][::step], Rss[::step, :], "c")
            ax.plot(ts[1:], Rtmedian, "r-", linewidth=2)
            ax.plot(ts[1:], RtLB, "r--", linewidth=2)
            ax.plot(ts[1:], RtUB, "r--", linewidth=2)
            ax.plot([0, ts[-1]], [1, 1], "k--", linewidth=2)
            ax.plot([timevect[-1], timevect[-1]], [0, 6], "k--", linewidth=2)
            ax.set_ylabel("$R_t$", fontsize=fs)
            ax.set_title("%s Ranked Model" % num2ordinal(rankx), fontsize=fs, fontweight="bold")
            ax.set_xlim(0, ts[-1])
            ax.set_ylim(0, 6)
            date_ticks(ax, datevecfirst1, timelags, DT, len(timevect2))

            ax = fig11.add_subplot(1, nModels, cc1)
            _interval_panel(ax, timevect2, LB1, UB1, np.median(curvesforecasts2, axis=1), data1, timevect[-1], fs)
            date_ticks(ax, datevecfirst1, timelags, DT, len(timevect2))
            ax.set_ylabel("%s %s" % (caddisease, datatype), fontsize=fs)
            ax.set_title("%s Ranked Model" % num2ordinal(rank1), fontsize=fs, fontweight="bold")

    if plot:
        for fig in (fig10, fig502, fig11):
            fig.set_facecolor("white")
            fig.tight_layout()

    # List the toolbox's warnings of this run once more.
    if verbose:
        subepidemicWarningLog("summary", "plotReproductionNumber")

    return np.array(param_Rt, dtype=float).reshape(-1, 5)


# <============================================================================>
# <=========================== Standalone curve plots =========================>
# <============================================================================>

_MODEL_NAMES = {
    0: "Generalized Growth Model",
    1: "Generalized Logistic Growth Model",
    2: "Generalized Richards Model",
    3: "Logistic Growth Model",
    4: "Richards Model",
    5: "Gompertz Model",
}


def plot_nsubepidemic(
    flag1_pass, r_pass, p_pass, a_pass, K_pass, npatches_pass, onset_thr_pass, I0_pass, windowsize1_pass,
    onset_fixed_pass=None, plot=True, ax=None,
):  # fmt: skip
    """Simulate and plot an n-sub-epidemic curve for given parameter values.

    Port of ``plot_nsubepidemic.m``.

    ``flag1_pass``        growth model, 0 to 5
    ``r_pass``, ``p_pass``, ``a_pass``, ``K_pass``  one value per sub-epidemic;
                          None (or empty) for a parameter the model does not
                          have (it is set to 1)
    ``npatches_pass``     number of sub-epidemics
    ``onset_thr_pass``    onset threshold C_thr
    ``I0_pass``           initial value C(0)
    ``windowsize1_pass``  duration; the time grid is 0, 1, ..., duration
    ``onset_fixed_pass``  0 asynchronous, 1 synchronous. When omitted, the
                          historical rule of the standalone function applies:
                          synchronous if the threshold is 0 or less.

    Returns ``totinc``, the total incidence curve. The components are drawn in
    color and their sum in black.

    The MATLAB function draws the generalized logistic model only: for the
    other models it asks for an argument ``q_pass`` that does not exist and
    returns without a curve. The port draws every model.
    """

    def given(v):
        return v is not None and np.size(v) > 0

    if not given(flag1_pass):
        raise SubEpiPredictError("SubEpiPredict:PlotInput", "flag1 value is required")
    flag1 = int(flag1_pass)
    if flag1 not in _MODEL_NAMES:
        raise SubEpiPredictError("SubEpiPredict:InvalidModelFlag", "flag1 must be an integer from 0 to 5.")
    if not given(npatches_pass):
        raise SubEpiPredictError("SubEpiPredict:PlotInput", "npatches value is required")
    npatches = int(npatches_pass)
    if not given(onset_thr_pass):
        raise SubEpiPredictError("SubEpiPredict:PlotInput", "onset_thr value is required")
    onset_thr = float(onset_thr_pass)

    # <==================================================================================>
    # <========================== Parameter values ======================================>
    # <==================================================================================>
    needs = {0: "rp", 1: "rpK", 2: "rpaK", 3: "rK", 4: "raK", 5: "rK"}[flag1]
    values = {"r": r_pass, "p": p_pass, "a": a_pass, "K": K_pass}
    params = {}
    for letter in "rpaK":
        if letter in needs:
            if not given(values[letter]):
                raise SubEpiPredictError("SubEpiPredict:PlotInput", "%s value is required" % letter)
            v = np.asarray(values[letter], dtype=float).ravel()
            if v.size == 1:
                v = np.repeat(v, npatches)
            params[letter] = v
        else:
            params[letter] = np.ones(npatches)

    # <==================================================================================>
    # <========================== Initial condition C(0) ================================>
    # <==================================================================================>
    if not given(I0_pass):
        raise SubEpiPredictError("SubEpiPredict:PlotInput", "Initial value C(0) is required")
    I0 = float(I0_pass)

    # <==================================================================================>
    # <========================== Time span =============================================>
    # <==================================================================================>
    if not given(windowsize1_pass):
        raise SubEpiPredictError("SubEpiPredict:PlotInput", "Simulation duration is required")
    timevect = np.arange(0, int(windowsize1_pass) + 1, dtype=float)

    if onset_fixed_pass is not None and np.size(onset_fixed_pass) > 0:
        onsetMode = int(onset_fixed_pass)
    else:
        # Legacy standalone interface only; fitted workflows always pass the mode.
        onsetMode = int(onset_thr <= 0)

    _, x, totinc, patchIncidence, _ = simulateSubepidemic(
        timevect, I0, npatches, onsetMode, onset_thr, flag1, params["r"], params["p"], params["a"], params["K"]
    )

    if plot:
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()
        fs = GetAdjustedFontSize()
        for j in range(npatches):
            ax.plot(timevect, patchIncidence[:, j], COMPONENT_COLORS[j % 7] + "-", linewidth=2)
        ax.plot(timevect, totinc, "k", linewidth=3)
        ax.set_ylabel("dC(t)/dt", fontsize=fs)
        ax.set_xlabel("Time", fontsize=fs)
        ax.set_xlim(0, timevect[-1])
        ax.set_ylim(0, float(np.max(totinc)) + 10)
        ax.figure.set_facecolor("white")

    return totinc


def plot_nsubepidemic_profiles(plot=True):
    """Four profiles of a three-sub-epidemic curve with different onset thresholds.

    Port of the script ``plot_nsubepidemic_profiles.m``: generalized logistic
    sub-epidemics with r = 0.18, p = 0.98 and K = 10000, 5000, 1000, C(0) = 1,
    220 time units, and C_thr = 100, 400, 3000, 5000 (panels A to D of
    figure 1). Returns the four total incidence curves as columns.
    """
    legends = ["A", "B", "C", "D"]
    onset_thr = [100, 400, 3000, 5000]
    duration1 = 220
    npatches = 3
    flag1 = 1
    rs = [0.18, 0.18, 0.18]
    ps = [0.98, 0.98, 0.98]
    as_ = None
    Ks = [10000, 5000, 1000]
    I0 = 1

    curves = []
    fig = None
    if plot:
        import matplotlib.pyplot as plt

        plt.close("all")
        fig = plt.figure(1, figsize=(13, 9))
    for cc1 in range(4):
        ax = fig.add_subplot(2, 2, cc1 + 1) if plot else None
        curve = plot_nsubepidemic(flag1, rs, ps, as_, Ks, npatches, onset_thr[cc1], I0, duration1, plot=plot, ax=ax)
        curves.append(curve)
        if plot:
            ax.set_title("$C_{thr}$=%s" % num2str(onset_thr[cc1]), fontsize=GetAdjustedFontSize(), fontweight="bold")
            ax.text(5, float(np.max(curve)) - 5, legends[cc1], fontsize=20, va="top")
    if plot:
        fig.tight_layout()
    return np.column_stack(curves)
