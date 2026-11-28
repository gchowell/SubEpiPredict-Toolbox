"""The fit driver: ``Run_Fit_subepidemicFramework``."""

from __future__ import annotations

import os

import numpy as np

from .aicc import get_nparams
from .errors import SubEpiPredictError, subepidemicWarning, subepidemicWarningLog
from .fitting import fittingModifiedLogisticFunctionPatchABC, fittingModifiedLogisticFunctionPatchMultiple
from .io import OutputNames, load_area_series, save_mat, to_json
from .noise import getMeanVarianceRatio
from .optionsets import options as default_options
from .optionsets import resolve_options
from .scores import subepidemicEvidenceRatios
from .settings import (
    resolve_fit_settings,
    subepidemicCheckSettings,
    subepidemicFitSignature,
    subepidemicSeed,
)
from .utils import Struct, cadtemporal_of, data_file_name, load_matrix, mat2str, resolve_input_file


def _figure_font():
    return 14


def Run_Fit_subepidemicFramework(
    outbreakx_pass=None,
    caddate1_pass=None,
    options=None,
    fit_settings=None,
    input_dir="./input",
    output_dir="./output",
    plot=True,
    verbose=True,
):
    """Fit and rank n-sub-epidemic models and bootstrap the top-ranked ones.

    Port of ``Run_Fit_subepidemicFramework.m``.

    Parameters
    ----------
    outbreakx_pass : area (column of the data file, from 1) or a list of
        areas; None = ``options.outbreakx``
    caddate1_pass : date stamp of the data file, ``'mm-dd-yyyy'``; None =
        ``options.caddate1``
    options : an :class:`~subepipredict.Options` object or a function
        returning one (default :func:`subepipredict.options`, the
        quick start). MATLAB reads ``options.m``.
    fit_settings : fit settings: a dict of overrides, a function returning
        one, or the result of :func:`subepidemicFitSettings`. MATLAB reads
        ``options_fit_settings.m``.
    input_dir : folder of the data files (default ``./input``; the example
        datasets of the repository are found from any folder)
    output_dir : folder of the result files (default ``./output``; created)
    plot : draw the figures of the fit (figures 10, 99, 100)
    verbose : print the progress and the summaries

    For every area the function loads the series, keeps the last
    ``calibrationperiod1`` points, drops the leading zeros, ranks the models
    over the number of sub-epidemics and the onset threshold by AICc
    (:func:`fittingModifiedLogisticFunctionPatchABC`) and bootstraps the
    ``topmodelsx`` best ones (:func:`fittingModifiedLogisticFunctionPatchMultiple`).

    Files written to ``output_dir`` (names as in MATLAB):
    ``ABC-ensem-...mat`` (the ranking) and
    ``modifiedLogisticPatch-ensem-...-rank-<k>.mat`` for each top model.

    Returns a list with one Struct per fitted area: ``outbreakx``, ``RMSES``,
    ``PS``, ``relativelik_i``, ``evidenceRatio``, ``topmodelsx``, ``screen``,
    ``ranked`` (list of Structs with the outputs of the bootstrap) and
    ``files``. (The MATLAB function returns nothing.)
    """
    opts = resolve_options(options, default_options)
    fitSettings = resolve_fit_settings(fit_settings)

    # <============================================================================>
    # <================================ Datasets properties =======================>
    # <============================================================================>
    outbreakx = opts.outbreakx if outbreakx_pass is None or np.size(outbreakx_pass) == 0 else outbreakx_pass
    caddate1 = opts.caddate1 if caddate1_pass is None or caddate1_pass == "" else caddate1_pass

    method1 = opts.method1
    dist1 = opts.dist1
    npatches_fixed = opts.npatches_fixed
    onset_fixed = opts.onset_fixed
    smoothfactor1 = opts.smoothfactor1
    calibrationperiod1 = opts.calibrationperiod1
    factor1 = None

    # <============================================================================>
    # <== Check the settings before anything else: every problem is reported at once ==>
    # <============================================================================>
    checkedSettings = opts.as_dict()
    checkedSettings["outbreakx"] = outbreakx
    subepidemicCheckSettings(checkedSettings)

    # The settings in one line, and a fresh list of warnings for this run.
    subepidemicWarningLog("start")
    if verbose:
        print(
            "Settings: method1 = %d, dist1 = %d, flag1 = %d, up to %d sub-epidemic(s), onset_fixed = %d, "
            "smoothfactor1 = %d, calibrationperiod1 = %d, %d start points, %d top models, B = %d."
            % (method1, dist1, opts.flag1, npatches_fixed, onset_fixed, smoothfactor1, calibrationperiod1,
               opts.numstartpoints, opts.topmodelsx, opts.B)
        )  # fmt: skip

    cadregion = opts.cadregion
    caddisease = opts.caddisease
    datatype = opts.datatype
    DT = opts.DT
    cumulative1 = opts.cumulative1
    cadtemporal = cadtemporal_of(DT)

    # Name of the file containing the time series data (rows=time, cols=regions)
    datafilename1 = data_file_name(cumulative1, cadtemporal, caddisease, datatype, cadregion, caddate1)

    numstartpoints = opts.numstartpoints
    topmodelsx = opts.topmodelsx
    M = opts.B  # number of bootstrap realizations to characterize parameter uncertainty

    # <==============================================================================>
    # <========================= Mathematical model =================================>
    # <==============================================================================>
    if npatches_fixed == 1:
        topmodelsx = 1
    flag1 = opts.flag1

    # <===========================================================================================================>
    # <====== Check that the number of estimated parameters is smaller than the number of data points ============>
    # <===========================================================================================================>
    numparams = get_nparams(method1, dist1, npatches_fixed, flag1, 1, onset_fixed)
    if numparams >= calibrationperiod1:
        raise SubEpiPredictError(
            "SubEpiPredict:TooManyParameters",
            "Number of estimated parameters should be smaller than the calibration period. Consider reducing "
            "<npatches_fixed> or increasing the length of the calibration period.",
        )

    # <==============================================================================>
    # <============ Load data and proceed to parameter estimation ===================>
    # <==============================================================================>
    datapath = resolve_input_file(datafilename1, input_dir)
    data = load_matrix(datapath)  # time series dataset (rows=time, cols=regions/areas)
    if data.size == 0:
        raise SubEpiPredictError("SubEpiPredict:EmptyData", "The dataset is empty")

    # Several areas: each area is fitted and saved as if the driver had been
    # called for that area alone.
    outbreakList = [int(v) for v in np.asarray(outbreakx).ravel()]
    areasSkipped = []  # areas with too few data points
    topmodelsRequested = topmodelsx  # the setting; topmodelsx can be lowered for one area

    # Every requested area must be a column of the data file; checked for all
    # areas now, not after the earlier ones have been fitted.
    if any(v > data.shape[1] for v in outbreakList):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidSettings",
            "outbreakx = %s, but the data file %s has %d column(s)."
            % (mat2str(outbreakList), datafilename1, data.shape[1]),
        )

    os.makedirs(output_dir, exist_ok=True)
    results = []

    for outbreak1 in outbreakList:

        if plot:
            import matplotlib.pyplot as plt

            plt.close("all")

        topmodelsx = topmodelsRequested

        # Random seed (setting rngSeed) and generator state, per area, so that an
        # area fitted in a call with several areas equals the same area fitted
        # alone. rngInfo is saved with the results of each ranked model.
        rngInfo = subepidemicSeed(
            "Run_Fit_subepidemicFramework, area %d" % outbreak1, fitSettings, verbose=verbose
        )

        names = OutputNames(opts, outbreak1, caddate1, output_dir)
        datafilename2 = names.cadfilename2 + ".txt"

        data1, _ = load_area_series(datafilename1, outbreak1, cumulative1, input_dir)

        # <==================================================================================>
        # <============ select length of the calibration period =============================>
        # <==================================================================================>
        if calibrationperiod1 < data1.size:
            data1[: data1.size - calibrationperiod1] = 0

        # <=============================================================================================>
        # <========== remove early zeros until a non-zero data point is found in the time series =======>
        # <=============================================================================================>
        index1 = np.flatnonzero(data1 > 0)
        if index1.size == 0:
            raise SubEpiPredictError(
                "SubEpiPredict:NoPositiveData",
                "Area %d has no positive value in the calibration period; nothing can be fitted." % outbreak1,
            )
        data1 = data1[index1[0] :]
        timelags = int(index1[0])  # keeps track of the epidemic onset timing (nonzero case)

        # <=============================================================================================>
        # <================ If time series is too short, do not proceed further ========================>
        # <=============================================================================================>
        if data1.size < 10 and npatches_fixed > 1:
            subepidemicWarning(
                "SubEpiPredict:AreaSkipped",
                "Area %d has %d data points after the calibration window and the leading zeros are removed. "
                "At least 10 are needed with npatches_fixed > 1. The area was skipped: nothing was fitted "
                "or saved for it.",
                outbreak1, data1.size,
            )  # fmt: skip
            areasSkipped.append(outbreak1)
            continue

        if method1 == 0 and dist1 == 2:
            # estimate <factor1>, the variance to mean ratio to scale the NB error structure
            binsize1 = 7
            ratios, _ = getMeanVarianceRatio(data1, binsize1, 2)
            keep = ratios[:, 0] > 0
            factor1 = float(np.mean(ratios[keep, 0])) if np.any(keep) else float("nan")

        # <==============================================================================================>
        # <================ Adjusted incidence data for further analysis ================================>
        # <==============================================================================================>
        data1 = np.column_stack([np.arange(data1.size, dtype=float), data1])
        data = data1.copy()
        t_window = np.arange(1, data.shape[0] + 1)
        data1 = data[t_window - 1, :]

        if plot:
            fig = plt.figure(1)
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(data1[:, 0], data1[:, 1], "ko", markerfacecolor="none")

        # <=======================================================================================>
        # <===================== Derive the best fitting sub-epidemic models =====================>
        # <=======================================================================================>
        RMSESx, PS, npatches, onset_thr, _, screen = fittingModifiedLogisticFunctionPatchABC(
            datafilename2, data1, DT, t_window, M, flag1, numstartpoints,
            method1=method1, dist1=dist1, npatches_fixed=npatches_fixed, onset_fixed=onset_fixed,
            smoothfactor1=smoothfactor1, calibrationperiod1=calibrationperiod1, fit_settings=fitSettings,
            output_dir=output_dir, plot=plot, verbose=verbose,
        )  # fmt: skip
        if screen.factor1 is not None:
            factor1 = screen.factor1

        # The screen can retain fewer models than topmodelsx (few onset
        # thresholds, candidates too large for the series, rejected fits).
        if topmodelsx > RMSESx.shape[0]:
            subepidemicWarning(
                "SubEpiPredict:TopModelsCapped",
                "The options ask for %d top-ranked models, but the screen of area %d retained %d. "
                "Only those %d are bootstrapped.",
                topmodelsx, outbreak1, RMSESx.shape[0], RMSESx.shape[0],
            )  # fmt: skip
            topmodelsx = RMSESx.shape[0]

        AICmin = RMSESx[0, 2]
        relativelik_i = np.exp((AICmin - RMSESx[:, 2]) / 2)

        # Evidence ratios of the top models against ranked model 1, with a
        # warning for those the data hardly support (report only).
        # What made this fit (revision and settings); saved with every ranked
        # file and compared by subepidemicRollingForecast before a fit is reused.
        fitSignature = subepidemicFitSignature(settings=fitSettings)
        evidenceRatio, _ = subepidemicEvidenceRatios(
            RMSESx[:, 2], topmodelsx, outbreak1, fitSettings.evidenceRatioLimit
        )
        if verbose:
            print(
                "Evidence ratios of ranked models 1 to %d against ranked model 1: %s"
                % (topmodelsx, "  ".join("%.3g" % v for v in evidenceRatio))
            )

        if plot:  # plot relative likelihoods of the models
            fig = plt.figure(99, figsize=(15, 4.5))
            fig.clf()
            ranks = np.arange(1, topmodelsx + 1)
            for k, (values, label) in enumerate(
                ((RMSESx[:topmodelsx, 2], "AICc"), (relativelik_i[:topmodelsx], "Relative likelihood"),
                 (evidenceRatio, "Evidence ratio"))
            ):  # fmt: skip
                ax = fig.add_subplot(1, 3, k + 1)
                ax.plot(ranks, values, "ko-", linewidth=2)
                ax.set_xlabel(r"$i_{th}$ ranked model", fontsize=_figure_font())
                ax.set_ylabel(label, fontsize=_figure_font())
                ax.tick_params(labelsize=_figure_font())
            fig.set_facecolor("white")
            fig.tight_layout()

        # <===========================================================================================>
        # <=============== Derive uncertainty for the <topmodelsx> best fitting models ===============>
        # <===========================================================================================>
        ranked = []
        files = [screen.file]
        for rank1 in range(1, topmodelsx + 1):
            (Phatss, npatches, onset_thr, curves, bestfit, data1b, P0, AICc_best, factor1_rank, d,
             bootstrapDiagnostics) = fittingModifiedLogisticFunctionPatchMultiple(
                RMSESx, PS, data1, DT, t_window, M, flag1, numstartpoints, rank1,
                method1=method1, dist1=dist1, onset_fixed=onset_fixed, smoothfactor1=smoothfactor1,
                factor1=factor1, fit_settings=fitSettings, plot=plot, verbose=verbose,
            )  # fmt: skip

            cadfilename1 = names.ranked_mat(rank1)
            save_mat(
                cadfilename1,
                dict(
                    Phatss=Phatss, npatches=npatches, onset_thr=onset_thr, curves=curves,
                    bestfit=bestfit.reshape(-1, 1), data1=data1b, data=data, P0=P0, AICc_best=AICc_best,
                    factor1=factor1_rank, d=d, bootstrapDiagnostics=bootstrapDiagnostics,
                    RMSESx=RMSESx, PS=PS, relativelik_i=relativelik_i.reshape(-1, 1),
                    evidenceRatio=evidenceRatio.reshape(-1, 1), t_window=t_window, timelags=timelags,
                    rank1=rank1, flag1=flag1, method1=method1, dist1=dist1, npatches_fixed=npatches_fixed,
                    onset_fixed=onset_fixed, smoothfactor1=smoothfactor1, calibrationperiod1=calibrationperiod1,
                    M=M, numstartpoints=numstartpoints, topmodelsx=topmodelsx,
                    topmodelsRequested=topmodelsRequested, outbreak1=outbreak1, outbreakx=outbreak1,
                    caddate1=caddate1, cadregion=cadregion, caddisease=caddisease, datatype=datatype, DT=DT,
                    cadtemporal=cadtemporal, cumulative1=cumulative1, datevecfirst1=opts.datevecfirst1,
                    datevecend1=opts.datevecend1, datafilename1=datafilename1, datafilename2=datafilename2,
                    fitSignature=fitSignature, fitSignatureJSON=to_json(fitSignature),
                    rngInfo=dict(seed=rngInfo.seed), madeBy="SubEpiPredict Python port",
                ),
            )  # fmt: skip
            files.append(cadfilename1)
            ranked.append(
                Struct(
                    rank=rank1, Phatss=Phatss, npatches=npatches, onset_thr=onset_thr, curves=curves,
                    bestfit=bestfit, data1=data1b, P0=P0, AICc_best=AICc_best, factor1=factor1_rank, d=d,
                    bootstrapDiagnostics=bootstrapDiagnostics, file=cadfilename1,
                )  # fmt: skip
            )
            if verbose:
                print("Ranked model %d complete." % rank1)

        results.append(
            Struct(
                outbreakx=outbreak1, caddate1=caddate1, RMSES=RMSESx, PS=PS, relativelik_i=relativelik_i,
                evidenceRatio=evidenceRatio, topmodelsx=topmodelsx, timelags=timelags, data1=data1,
                screen=screen, ranked=ranked, files=files, rngInfo=rngInfo, fitSignature=fitSignature,
            )  # fmt: skip
        )

    if areasSkipped and verbose:
        print(
            "Areas skipped for too few data points (nothing was saved for them): %s" % mat2str(areasSkipped)
        )

    # List the toolbox's warnings of this run once more.
    if verbose:
        subepidemicWarningLog("summary", "Run_Fit_subepidemicFramework")

    return results
