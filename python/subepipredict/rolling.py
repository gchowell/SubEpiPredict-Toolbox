"""Rolling-origin evaluation: fit and forecast from a list of dates and score
every horizon.

Ports of ``subepidemicRollingForecast.m`` and ``subepidemicRollingSummary.m``.
"""

from __future__ import annotations

import contextlib
import glob
import io
import math
import os
import shutil
import time
import warnings

import numpy as np

from .errors import SubEpiPredictError, SubEpiPredictWarning, subepidemicWarningLog
from .forecast import plotForecast_subepidemicFramework
from .io import OutputNames, from_json, load_mat, scalar, writetable
from .optionsets import options as default_options
from .optionsets import options_forecast as default_options_forecast
from .optionsets import resolve_options
from .runners import Run_Fit_subepidemicFramework
from .settings import resolve_fit_settings, subepidemicFitSignature
from .utils import Struct, cadtemporal_of, data_file_name, datenum, datestr, datevec, find_input_file, load_matrix, mat2str

SUMMARY_COLUMNS = [
    "model", "horizon", "dates", "MAE", "Coverage 95%PI", "WIS", "WIS relative to baseline (ratio of means)",
    "WIS relative to baseline (typical ratio)", "share of dates better than baseline",
    "MAE relative to baseline (ratio of means)", "WIS at this step", "WIS at this step relative to baseline",
    "Coverage 95%PI at this step",
]  # fmt: skip

SCORES_COLUMNS = ["date number", "model", "horizon", "MAE", "MSE", "Coverage 95%PI", "WIS", "WIS of the baseline"]

STATUS_COLUMNS = ["date number", "data rows", "input file (1 written, 2 existing, 3 existing and different)",
                  "fit (1 new, 2 reused)", "scored", "kinds of warning", "seconds"]  # fmt: skip


def subepidemicRollingSummary(scores, horizon):
    """Average forecast scores over dates, by model and horizon.

    Port of ``subepidemicRollingSummary.m``.

    ``scores``   one row per date, model and horizon, as collected by
                 :func:`subepidemicRollingForecast` (``R.scores``): date
                 number, model, horizon, MAE, MSE, coverage of the 95% PI (%),
                 WIS, WIS of the baseline. Model 0 is the baseline. The scores
                 of horizon h are averages over the first h time units.
    ``horizon``  the largest horizon

    Returns ``(summary, columns)``. One row per model and horizon; columns
    (numbered from 1): 1 model; 2 horizon h; 3 number of dates with a score;
    4 MAE over the first h time units, mean over the dates; 5 coverage of the
    95% PI over the first h time units (%); 6 WIS over the first h time units;
    7 WIS relative to the baseline, ratio of means; 8 the same, typical ratio
    (geometric mean over the dates of the ratio); 9 share of dates on which
    the model's WIS is below the baseline's (a tie counts half); 10 MAE
    relative to the baseline, ratio of means; 11 WIS at step h alone; 12 WIS
    at step h alone relative to the baseline; 13 coverage at step h alone.

    The score at step h alone is recovered from the averages:
    ``h*score(h) - (h-1)*score(h-1)``.
    """
    columns = list(SUMMARY_COLUMNS)
    summary = np.zeros((0, 13))
    scores = np.asarray(scores, dtype=float)
    if scores.size == 0:
        return summary, columns
    scores = scores.reshape(-1, 8)
    horizon = int(horizon)
    dates = np.unique(scores[:, 0])
    models = np.unique(scores[:, 1])
    nDates = dates.size

    # scores of each model as date-by-horizon arrays (NaN where there is no row)
    maeOf, covOf, wisOf = [], [], []
    for model in models:
        mae = np.full((nDates, horizon), np.nan)
        cov = np.full((nDates, horizon), np.nan)
        wis = np.full((nDates, horizon), np.nan)
        rowsModel = scores[scores[:, 1] == model, :]
        dateIndex = np.searchsorted(dates, rowsModel[:, 0])
        for k in range(rowsModel.shape[0]):
            h = int(rowsModel[k, 2])
            if 1 <= h <= horizon:
                mae[dateIndex[k], h - 1] = rowsModel[k, 3]
                cov[dateIndex[k], h - 1] = rowsModel[k, 5]
                wis[dateIndex[k], h - 1] = rowsModel[k, 6]
        maeOf.append(mae)
        covOf.append(cov)
        wisOf.append(wis)

    baseIdx = np.flatnonzero(models == 0)
    base = int(baseIdx[0]) if baseIdx.size else None
    steps = np.tile(np.arange(1, horizon + 1, dtype=float), (nDates, 1))

    def atStep(avg):
        return steps * avg - (steps - 1) * np.column_stack([np.zeros(nDates), avg[:, :-1]])

    rows = []
    for m, model in enumerate(models):
        maeM, covM, wisM = maeOf[m], covOf[m], wisOf[m]
        # scores at step h alone
        wisStep = atStep(wisM)
        covStep = atStep(covM)
        if base is None:
            maeB = np.full((nDates, horizon), np.nan)
            wisB = maeB
            wisStepB = maeB
        else:
            maeB = maeOf[base]
            wisB = wisOf[base]
            wisStepB = atStep(wisB)
        for h in range(horizon):
            row = np.full(13, np.nan)
            row[0] = model
            row[1] = h + 1
            have = np.isfinite(maeM[:, h])
            row[2] = np.sum(have)
            if row[2] == 0:
                continue
            row[3] = np.mean(maeM[have, h])
            haveCov = np.isfinite(covM[:, h])
            if np.any(haveCov):
                row[4] = np.mean(covM[haveCov, h])
            haveWIS = np.isfinite(wisM[:, h])
            if np.any(haveWIS):
                row[5] = np.mean(wisM[haveWIS, h])
            both = haveWIS & np.isfinite(wisB[:, h])
            if np.any(both):
                if np.sum(wisB[both, h]) > 0:
                    row[6] = np.sum(wisM[both, h]) / np.sum(wisB[both, h])
                positive = both & (wisM[:, h] > 0) & (wisB[:, h] > 0)
                if np.any(positive):
                    row[7] = math.exp(np.mean(np.log(wisM[positive, h] / wisB[positive, h])))
                if model != 0:
                    row[8] = (
                        np.sum(wisM[both, h] < wisB[both, h]) + 0.5 * np.sum(wisM[both, h] == wisB[both, h])
                    ) / np.sum(both)
            bothMAE = have & np.isfinite(maeB[:, h])
            if np.any(bothMAE) and np.sum(maeB[bothMAE, h]) > 0:
                row[9] = np.sum(maeM[bothMAE, h]) / np.sum(maeB[bothMAE, h])
            haveStep = np.isfinite(wisStep[:, h])
            if np.any(haveStep):
                row[10] = np.mean(wisStep[haveStep, h])
            bothStep = haveStep & np.isfinite(wisStepB[:, h])
            if np.any(bothStep) and np.sum(wisStepB[bothStep, h]) > 0:
                row[11] = np.sum(wisStep[bothStep, h]) / np.sum(wisStepB[bothStep, h])
            haveCovStep = np.isfinite(covStep[:, h])
            if np.any(haveCovStep):
                row[12] = np.mean(covStep[haveCovStep, h])
            rows.append(row)
    if rows:
        summary = np.vstack(rows)
    return summary, columns


def _parse_dates(forecastDates):
    """Date numbers of the forecast dates, in increasing order, without repeats."""
    if isinstance(forecastDates, str):
        forecastDates = [forecastDates]
    forecastDates = list(np.asarray(forecastDates, dtype=object).ravel()) if not isinstance(forecastDates, list) else forecastDates
    dateNumbers = []
    for k, item in enumerate(forecastDates, start=1):
        if isinstance(item, str):
            d = item.strip()
            digits = d[0:2] + d[3:5] + d[6:10] if len(d) == 10 else ""
            if len(d) != 10 or d[2] != "-" or d[5] != "-" or not digits.isdigit():
                raise SubEpiPredictError(
                    "SubEpiPredict:RollingInput",
                    "Forecast date %d is '%s'. Dates are written mm-dd-yyyy, e.g. 05-11-2020." % (k, d),
                )
            month, day, year = int(d[0:2]), int(d[3:5]), int(d[6:10])
            try:
                import datetime

                datetime.date(year, month, day)
            except ValueError:
                raise SubEpiPredictError(
                    "SubEpiPredict:RollingInput",
                    "Forecast date %d, '%s', is not a calendar date (dates are written mm-dd-yyyy)." % (k, d),
                ) from None
            dateNumbers.append(datenum(year, month, day))
        else:
            dateNumbers.append(int(round(float(item))))
    return np.unique(np.array(dateNumbers, dtype=int))


def subepidemicRollingForecast(
    forecastDates,
    forecastingperiod=None,
    outbreakx=None,
    weight_type1=None,
    *args,
    refit=False,
    closeFigures=True,
    quiet=False,
    keepReplaced=True,
    options=None,
    options_forecast=None,
    fit_settings=None,
    input_dir="./input",
    output_dir="./output",
    plot=True,
):
    """Fit and forecast from a list of dates and score every horizon.

    Port of ``subepidemicRollingForecast.m``. The function repeats the fit
    (:func:`Run_Fit_subepidemicFramework`) and the scored forecast
    (:func:`plotForecast_subepidemicFramework`) from every date of a list and
    averages the scores of each horizon over the dates, next to the flat
    baseline.

    Parameters
    ----------
    forecastDates : list of dates ``'mm-dd-yyyy'`` (the last date of the data
        used for each fit), or of MATLAB serial date numbers
        (:func:`subepipredict.utils.datenum`)
    forecastingperiod : horizon in time units (None = ``options_forecast``)
    outbreakx : one area, a column of the data file (None = the options)
    weight_type1 : ensemble weights (None = ``options_forecast``)
    refit : False (default): a date whose fit is already in ``output_dir``,
        made with the same B, numstartpoints and topmodelsx, by the same fit
        revision with the same fit settings and complete, is not fitted again.
        True: every date is fitted again.
    closeFigures : True (default): the figures a date opens are closed before
        the next date
    quiet : True: the text of each fit and forecast goes to a log file per
        date in ``output_dir`` (``rolling-log-...txt``) and each date prints a
        few lines
    keepReplaced : True (default): before a date is fitted again, the fit
        files of that date already in ``output_dir`` are copied to
        ``output_dir/replaced-by-rolling/``, unless a copy is there already
    options, options_forecast, fit_settings, input_dir, output_dir, plot :
        as in the drivers

    The MATLAB name-value pairs can also be given as in MATLAB:
    ``subepidemicRollingForecast(dates, 4, 52, -1, 'quiet', True)``.

    Data. The data of every date come from one file: the file of the most
    recent date (``datevecend1``), which also supplies the later observations
    the forecasts are scored against. For a date without its own file, the
    file is written to ``input_dir`` from the first rows of that most recent
    file. A file that already exists is used as it is, and a difference from
    the most recent file is reported.

    Returns a Struct ``R``: ``dates``, ``status`` (one row per date, see
    ``statusColumns``), ``notes``, ``scores`` (one row per date, model and
    horizon), ``summary`` (one row per model and horizon,
    :func:`subepidemicRollingSummary`), the column names, ``horizon`` and the
    names of the two files written (``scoresFile``, ``summaryFile``). Model
    codes: 0 baseline, k ranked k, 100+k Ensemble(k).
    """
    # ---------------------------------------------------------------- inputs
    if forecastDates is None or np.size(forecastDates) == 0:
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "Give the forecast dates: a list of 'mm-dd-yyyy' or of serial date numbers.",
        )
    if len(args) % 2 != 0:
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput", "Options are given as name and value pairs, e.g. 'refit',True."
        )
    for k in range(0, len(args), 2):
        name = str(args[k]).lower()
        if name == "refit":
            refit = bool(args[k + 1])
        elif name == "closefigures":
            closeFigures = bool(args[k + 1])
        elif name == "quiet":
            quiet = bool(args[k + 1])
        elif name == "keepreplaced":
            keepReplaced = bool(args[k + 1])
        else:
            raise SubEpiPredictError(
                "SubEpiPredict:RollingInput",
                "Unknown option '%s'. The options are 'refit', 'closeFigures', 'quiet' and 'keepReplaced'."
                % args[k],
            )

    opts = resolve_options(options, default_options)
    fopts = resolve_options(options_forecast, default_options_forecast)
    fitSettings = resolve_fit_settings(fit_settings)

    cumulative1 = opts.cumulative1
    DT = opts.DT
    datevecfirst1 = opts.datevecfirst1
    datevecend1 = opts.datevecend1
    numstartpoints = opts.numstartpoints
    topmodelsSetting = opts.topmodelsx
    B = opts.B
    if opts.npatches_fixed == 1:
        topmodelsSetting = 1  # as in the fit driver

    if forecastingperiod is None:
        forecastingperiod = fopts.forecastingperiod
    if outbreakx is None:
        outbreakx = opts.outbreakx
    if weight_type1 is None:
        weight_type1 = fopts.weight_type1
    if np.size(outbreakx) != 1:
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "Give one area (outbreakx is a single column of the data file); run the function once per area.",
        )
    outbreakx = np.asarray(outbreakx).ravel()[0]
    if not (np.ndim(forecastingperiod) == 0 and forecastingperiod >= 1 and forecastingperiod == round(forecastingperiod)):
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "The forecast horizon must be a whole number of at least 1 (it is %s)." % mat2str(forecastingperiod),
        )
    forecastingperiod = int(forecastingperiod)
    if fopts.getperformance != 1:
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "options_forecast has getperformance = %d. Set it to 1: without it no forecast is scored."
            % fopts.getperformance,
        )
    if DT not in (1, 7, 365):
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "DT = %s in the options; the toolbox reads daily (1), weekly (7) and yearly (365) data." % mat2str(DT),
        )
    cadtemporal = cadtemporal_of(DT)

    def fileName(dateString):
        return data_file_name(cumulative1, cadtemporal, opts.caddisease, opts.datatype, opts.cadregion, dateString)

    dateNumbers = _parse_dates(forecastDates)
    nDates = dateNumbers.size

    # ---------------------------------------------- the most recent data file
    latestName = fileName(datestr(datenum(datevecend1), "mm-dd-yyyy"))
    latestPath = find_input_file(latestName, input_dir)
    if latestPath is None:
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "The file %s was not found in %s. It is named by datevecend1 in the options and holds the whole "
            "series: the data of every forecast date and the later observations the forecasts are scored against."
            % (latestName, input_dir),
        )
    latest = load_matrix(latestPath)
    if outbreakx < 1 or outbreakx > latest.shape[1] or outbreakx != round(outbreakx):
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "outbreakx = %s, but the data file %s has %d column(s)." % (mat2str(outbreakx), latestName, latest.shape[1]),
        )
    outbreakx = int(outbreakx)
    column = outbreakx - 1
    nLatest = latest.shape[0]

    # number of data rows up to each date; later data needed: rows+1 ... rows+horizon
    firstNumber = datenum(datevecfirst1)
    if DT == 365:
        rowsAtDate = np.array([datevec(d)[0] - int(datevecfirst1[0]) + 1 for d in dateNumbers])
        onGrid = np.ones(nDates, dtype=bool)
        steps = rowsAtDate - 1.0
    else:
        steps = (dateNumbers - firstNumber) / DT
        onGrid = np.abs(steps - np.round(steps)) < 1e-9
        rowsAtDate = np.round(steps).astype(int) + 1
    if np.any(~onGrid):
        k = int(np.flatnonzero(~onGrid)[0])
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "Forecast date %s is not one of the dates of the series: the series starts on %s and has one value "
            "every %d days. The nearest dates of the series are %s and %s."
            % (datestr(dateNumbers[k]), datestr(firstNumber), DT,
               datestr(firstNumber + int(math.floor(steps[k])) * DT), datestr(firstNumber + int(math.ceil(steps[k])) * DT)),
        )  # fmt: skip
    if np.any(rowsAtDate < 1):
        k = int(np.flatnonzero(rowsAtDate < 1)[0])
        raise SubEpiPredictError(
            "SubEpiPredict:RollingInput",
            "Forecast date %s is before the first date of the series (%s)."
            % (datestr(dateNumbers[k]), datestr(firstNumber)),
        )

    # status columns: 1 date number, 2 data rows, 3 input file (1 written, 2 existing and equal to the
    # most recent file, 3 existing and different), 4 fit (1 new, 2 reused), 5 scored (1/0), 6 warnings, 7 seconds
    status = np.column_stack([dateNumbers, rowsAtDate, np.full((nDates, 2), np.nan), np.zeros((nDates, 2)),
                              np.full(nDates, np.nan)]).astype(float)  # fmt: skip
    notes = ["" for _ in range(nDates)]
    print(
        "\nRolling forecasts: %d date(s) from %s to %s, horizon %d, area %d, %s %s %s (%s)."
        % (nDates, datestr(dateNumbers[0]), datestr(dateNumbers[-1]), forecastingperiod, outbreakx, cadtemporal,
           opts.caddisease, opts.datatype, opts.cadregion)
    )  # fmt: skip

    os.makedirs(output_dir, exist_ok=True)
    scores = np.zeros((0, 8))

    for k in range(nDates):
        dateString = datestr(dateNumbers[k])
        timeStart = time.time()
        print("\n===== Rolling forecast %d of %d: data to %s =====" % (k + 1, nDates, dateString))

        # ---- later data
        if rowsAtDate[k] + forecastingperiod > nLatest:
            notes[k] = "skipped: the forecast needs data up to row %d and %s has %d rows" % (
                rowsAtDate[k] + forecastingperiod, latestName, nLatest,
            )  # fmt: skip
            warnings.warn(SubEpiPredictWarning(
                "SubEpiPredict:RollingDateSkipped", "Forecast date %s %s." % (dateString, notes[k])), stacklevel=2)  # fmt: skip
            continue

        # ---- input file of this date
        dateName = fileName(dateString)
        datePath = find_input_file(dateName, input_dir)
        if datePath is not None:
            existing = load_matrix(datePath)
            if existing.shape[0] != rowsAtDate[k] or existing.shape[1] != latest.shape[1]:
                notes[k] = (
                    "skipped: the file %s has %d row(s) and %d column(s); a series that starts on %s needs "
                    "%d row(s) on that date and %d column(s)"
                    % (dateName, existing.shape[0], existing.shape[1], datestr(firstNumber), rowsAtDate[k],
                       latest.shape[1])
                )  # fmt: skip
                warnings.warn(SubEpiPredictWarning(
                    "SubEpiPredict:RollingDateSkipped", "Forecast date %s %s." % (dateString, notes[k])), stacklevel=2)  # fmt: skip
                continue
            if np.array_equal(existing[:, column], latest[: rowsAtDate[k], column]):
                status[k, 2] = 2
            else:
                status[k, 2] = 3
                nDifferent = int(np.sum(existing[:, column] != latest[: rowsAtDate[k], column]))
                notes[k] = (
                    "the existing file of this date differs from the most recent file in %d of %d values "
                    "(the fit uses the existing file)" % (nDifferent, rowsAtDate[k])
                )
                print("Data: %s exists and %s." % (dateName, notes[k]))
        else:
            os.makedirs(input_dir, exist_ok=True)
            try:
                with open(os.path.join(input_dir, dateName), "w") as fid:
                    for row in latest[: rowsAtDate[k], :]:
                        fid.write("\t".join("%.15g" % v for v in row) + "\n")
            except OSError:
                raise SubEpiPredictError(
                    "SubEpiPredict:RollingInput",
                    "The file %s could not be written." % os.path.join(input_dir, dateName),
                ) from None
            status[k, 2] = 1
            print("Data: %s written to %s (the first %d rows of %s)." % (dateName, input_dir, rowsAtDate[k], latestName))

        figuresBefore = None
        if plot:
            import matplotlib.pyplot as plt

            figuresBefore = set(plt.get_fignums())

        try:
            # ---- fit, unless a fit of this date with the same B and start points is saved
            names = OutputNames(opts, outbreakx, dateString, output_dir)
            fitName = names.ranked_mat(1)
            reuse = False
            if not refit and os.path.isfile(fitName):
                # The fit is used again if it was made with the same B, start
                # points and number of top models, and if it is complete: the file
                # of its last rank exists (an interrupted fit leaves the first ranks only).
                saved = load_mat(fitName)
                sameSettings = (
                    "M" in saved and "numstartpoints" in saved
                    and scalar(saved.M) == B and scalar(saved.numstartpoints) == numstartpoints
                )  # fmt: skip
                if sameSettings and "topmodelsRequested" in saved:
                    sameSettings = scalar(saved.topmodelsRequested) == topmodelsSetting
                lastRank = int(scalar(saved.topmodelsx)) if "topmodelsx" in saved else topmodelsSetting
                complete = os.path.isfile(names.ranked_mat(lastRank))
                # The fit must also have been made by the current fit revision
                # with the same fit settings.
                savedSignature = None
                if "fitSignatureJSON" in saved:
                    try:
                        savedSignature = from_json(saved.fitSignatureJSON)
                    except (TypeError, ValueError):
                        savedSignature = None
                sameSignature, signatureDifferences = subepidemicFitSignature(
                    savedSignature, settings=fitSettings, compare=True
                )
                del saved
                if sameSettings and complete and sameSignature:
                    reuse = True
                elif sameSettings and complete:
                    print(
                        "Fit: the saved fit of this date is not reused (%s); the date is fitted again."
                        % ", ".join(signatureDifferences)
                    )
                elif not sameSettings:
                    print(
                        "Fit: a saved fit of this date was made with other values of B, numstartpoints or "
                        "topmodelsx; the date is fitted again."
                    )
                else:
                    print(
                        "Fit: the saved fit of this date is incomplete (the file of rank %d is missing); "
                        "the date is fitted again." % lastRank
                    )

            nWarnings = 0
            consoleText = ""
            if reuse:
                status[k, 3] = 2
                print("Fit: the saved fit of this date is used (give refit=True to fit every date again).")
            else:
                # ---- keep a copy of the fit files this fit replaces
                if keepReplaced:
                    existingFiles = glob.glob(glob.escape(fitName[:-5]) + "*.mat")
                    if os.path.isfile(names.abc_mat()):
                        existingFiles.append(names.abc_mat())
                    nCopied = 0
                    for path in existingFiles:
                        folder = os.path.join(output_dir, "replaced-by-rolling")
                        os.makedirs(folder, exist_ok=True)
                        copyName = os.path.join(folder, os.path.basename(path))
                        if not os.path.isfile(copyName):
                            shutil.copyfile(path, copyName)
                            nCopied += 1
                    if nCopied > 0:
                        print(
                            "Fit: %d fit file(s) of this date already in %s were copied to %s before the new fit."
                            % (nCopied, output_dir, os.path.join(output_dir, "replaced-by-rolling"))
                        )
                kwargs = dict(options=opts, fit_settings=fitSettings, input_dir=input_dir, output_dir=output_dir,
                              plot=plot)  # fmt: skip
                if quiet:
                    buffer = io.StringIO()
                    with contextlib.redirect_stdout(buffer), warnings.catch_warnings():
                        warnings.simplefilter("always", SubEpiPredictWarning)
                        with contextlib.redirect_stderr(buffer):
                            Run_Fit_subepidemicFramework(outbreakx, dateString, **kwargs)
                    consoleText = buffer.getvalue()
                else:
                    Run_Fit_subepidemicFramework(outbreakx, dateString, **kwargs)
                status[k, 3] = 1
                nWarnings += len(subepidemicWarningLog("get"))
                if quiet:
                    print("Fit: new (%d kind(s) of warning)." % nWarnings)

            # ---- forecast and scores
            kwargs = dict(options=opts, options_forecast=fopts, fit_settings=fitSettings, input_dir=input_dir,
                          output_dir=output_dir, plot=plot)  # fmt: skip
            logName = ""
            if quiet:
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer), warnings.catch_warnings():
                    warnings.simplefilter("always", SubEpiPredictWarning)
                    with contextlib.redirect_stderr(buffer):
                        _, _, _, byHorizon = plotForecast_subepidemicFramework(
                            outbreakx, dateString, forecastingperiod, weight_type1, **kwargs
                        )
                consoleText += buffer.getvalue()
                logName = os.path.join(
                    output_dir,
                    "rolling-log-%s-%s-%s-%s-area-%d-%s.txt"
                    % (cadtemporal, opts.caddisease, opts.datatype, opts.cadregion, outbreakx, dateString),
                )
                try:
                    with open(logName, "w") as fid:
                        fid.write(consoleText)
                except OSError:
                    pass
            else:
                _, _, _, byHorizon = plotForecast_subepidemicFramework(
                    outbreakx, dateString, forecastingperiod, weight_type1, **kwargs
                )
            nWarningsForecast = len(subepidemicWarningLog("get"))
            nWarnings += nWarningsForecast
            status[k, 5] = nWarnings
            if quiet:
                scoreText = "no scores"
                if byHorizon.size:
                    rowBase = byHorizon[(byHorizon[:, 0] == 0) & (byHorizon[:, 1] == forecastingperiod), :]
                    rowOne = byHorizon[(byHorizon[:, 0] == 1) & (byHorizon[:, 1] == forecastingperiod), :]
                    if rowBase.size and rowOne.size:
                        scoreText = "WIS over %d time unit(s): ranked 1 %.2f, baseline %.2f, ratio %.2f" % (
                            forecastingperiod, rowOne[0, 5], rowBase[0, 5], rowOne[0, 5] / rowBase[0, 5],
                        )  # fmt: skip
                print(
                    "Forecast: %s (%d kind(s) of warning). %.0f s for this date.\nLog: %s"
                    % (scoreText, nWarningsForecast, time.time() - timeStart, logName)
                )
            if byHorizon.size == 0:
                notes[k] = (notes[k] + " no scores were returned").strip()
            else:
                # WIS of the baseline at the same horizon, next to every row
                baseRows = byHorizon[byHorizon[:, 0] == 0, :]
                baseWIS = np.full(forecastingperiod, np.nan)
                baseWIS[baseRows[:, 1].astype(int) - 1] = baseRows[:, 5]
                scores = np.vstack([
                    scores,
                    np.column_stack([np.full(byHorizon.shape[0], float(dateNumbers[k])), byHorizon[:, :6],
                                     baseWIS[byHorizon[:, 1].astype(int) - 1]]),
                ])  # fmt: skip
                status[k, 4] = 1
        except Exception as problem:
            # one date that fails does not stop the others
            message = getattr(problem, "message", None) or str(problem)
            notes[k] = "failed: %s" % message
            warnings.warn(SubEpiPredictWarning(
                "SubEpiPredict:RollingDateFailed",
                "Forecast date %s failed and is left out: %s" % (dateString, message)), stacklevel=2)  # fmt: skip

        if closeFigures and plot:
            import matplotlib.pyplot as plt

            for number in plt.get_fignums():
                if number not in figuresBefore:
                    plt.close(number)
        status[k, 6] = time.time() - timeStart

    # ----------------------------------------------------------- summary
    summary, summaryColumns = subepidemicRollingSummary(scores, forecastingperiod)

    R = Struct()
    R.dates = [datestr(d) for d in dateNumbers]
    R.status = status
    R.statusColumns = list(STATUS_COLUMNS)
    R.notes = notes
    R.scores = scores
    R.scoresColumns = list(SCORES_COLUMNS)
    R.summary = summary
    R.summaryColumns = summaryColumns
    R.horizon = forecastingperiod
    R.scoresFile = ""
    R.summaryFile = ""

    nScored = int(np.sum(status[:, 4] == 1))
    if nScored == 0:
        print("\nRolling forecasts: no date was scored.")
        for k in range(nDates):
            print("  %s  %s" % (R.dates[k], notes[k]))
        return R

    # files, named by the first and the last scored date
    fileTag = "onsetfixed-%d-flag1-%d-method-%d-dist-%d-horizon-%d-weight_type-%d-%s-%s-%s-%s-area-%d-%s-to-%s.csv" % (
        opts.onset_fixed, opts.flag1, opts.method1, opts.dist1, forecastingperiod, weight_type1, cadtemporal,
        opts.caddisease, opts.datatype, opts.cadregion, outbreakx, datestr(np.min(scores[:, 0])),
        datestr(np.max(scores[:, 0])),
    )  # fmt: skip
    dv = np.array([datevec(d) for d in scores[:, 0]], dtype=float)
    R.scoresFile = os.path.join(output_dir, "rolling-forecast-scores-" + fileTag)
    writetable(R.scoresFile, np.column_stack([dv, scores[:, 1:8]]),
               ["year", "month", "day", "model", "horizon", "MAE", "MSE", "Coverage 95%PI", "WIS",
                "WIS of the baseline"])  # fmt: skip
    R.summaryFile = os.path.join(output_dir, "rolling-forecast-summary-" + fileTag)
    writetable(R.summaryFile, summary, summaryColumns)

    # table on the screen
    if forecastingperiod <= 6:
        shown = list(range(1, forecastingperiod + 1))
    elif DT == 1 and forecastingperiod >= 14:
        shown = sorted(set(list(range(7, forecastingperiod + 1, 7)) + [forecastingperiod]))
    else:
        shown = sorted(set(int(round(v)) for v in np.linspace(1, forecastingperiod, 5)))
    models = np.unique(summary[:, 0])
    print(
        "\nRolling forecasts: %d of %d date(s) scored, %s to %s, horizon %d."
        % (nScored, nDates, R.dates[0], R.dates[-1], forecastingperiod)
    )
    blocks = [
        ("WIS over the first h time units, relative to the baseline (ratio of means; below 1 = better than the baseline)", 7),
        ("The same, typical ratio (geometric mean over the dates)", 8),
        ("Share of dates with a WIS over the first h time units below the baseline's", 9),
        ("WIS at step h alone, relative to the baseline (ratio of means)", 12),
        ("Coverage of the 95% prediction interval over the first h time units (%; nominal 95)", 5),
    ]
    for title, col in blocks:
        print("\n%s" % title)
        print("  %-14s %6s" % (" ", "dates") + "".join(" %8s" % ("h=%d" % h) for h in shown))
        for m in models:
            if m == 0 and col in (7, 8, 9, 12):
                continue
            if m == 0:
                modelLabel = "baseline"
            elif m > 100:
                modelLabel = "Ensemble(%d)" % (m - 100)
            else:
                modelLabel = "ranked %d" % m
            rowsModel = summary[summary[:, 0] == m, :]
            line = "  %-14s %6d" % (modelLabel, int(np.max(rowsModel[:, 2])))
            for h in shown:
                value = rowsModel[rowsModel[:, 1] == h, col - 1]
                value = float(value[0]) if value.size else float("nan")
                line += (" %8.1f" if col == 5 else " %8.2f") % value
            print(line)
    print("\nFiles: %s\n       %s" % (R.scoresFile, R.summaryFile))
    if np.any(status[:, 4] != 1) or np.any(status[:, 2] == 3):
        print("Notes:")
        for k in range(nDates):
            if notes[k]:
                print("  %s  %s" % (R.dates[k], notes[k]))
    if np.any(status[:, 5] > 0):
        print(
            "Dates with toolbox warnings (see the list printed after each run above): %s"
            % ", ".join(R.dates[k] for k in range(nDates) if status[k, 5] > 0)
        )
    scoredDates = dateNumbers[status[:, 4] == 1]
    if DT == 1 and nScored > 1 and np.all(scoredDates % 7 == scoredDates[0] % 7):
        print(
            "All scored dates fall on the same day of the week: with a weekly reporting cycle in the data, the "
            "scores at a single step h mix the horizon with the day of the week. Compare steps 7 days apart, or "
            "use the scores over the first 7, 14, ... days."
        )
    if nScored < 10:
        print(
            "Only %d date(s): differences between the models and the baseline are weak evidence at this number."
            % nScored
        )
    return R
