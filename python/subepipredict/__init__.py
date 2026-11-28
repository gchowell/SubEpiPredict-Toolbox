"""SubEpiPredict: ensemble n-sub-epidemic modeling and forecasting.

Python port of the MATLAB toolbox by Gerardo Chowell
(https://github.com/gchowell/SubEpiPredict-Toolbox). The functions keep their
MATLAB names and argument order::

    import subepipredict as sp

    sp.Run_Fit_subepidemicFramework(52, '05-11-2020')
    ranked, ensembles, baseline, byHorizon = sp.plotForecast_subepidemicFramework(52, '05-11-2020', 4, -1)

See the README of the ``python`` folder for the differences from MATLAB.
"""

from __future__ import annotations

from .aicc import get_nparams, getAICc
from .doubling import getDoublingTimeCurve
from .errors import SubEpiPredictError, SubEpiPredictWarning, subepidemicWarning, subepidemicWarningLog
from .figures import (
    GetAdjustedFontSize,
    subepidemicEstimateLine,
    subepidemicPanelGrid,
    subepidemicParameterLabel,
)
from .fitting import (
    fitSubepidemicModel,
    fittingModifiedLogisticFunctionPatchABC,
    fittingModifiedLogisticFunctionPatchMultiple,
    onsetThresholdGrid,
    shutdown_pool,
)
from .forecast import (
    getensemblesubepidemicmodels,
    plotForecast_subepidemicFramework,
    subepidemicAvailableTopModels,
    subepidemicSkippedRefits,
)
from .identifiability import (
    sci_interval,
    subepidemicBoundHitFractions,
    subepidemicBoundSummary,
    subepidemicIdentifiability,
    subepidemicPII,
    subepidemicPrintPII,
    subepidemicQuantileMCSE,
    subepidemicWritePII,
)
from .io import OutputNames, getData, load_mat, readtable, save_mat, writetable
from .likelihood import subepidemicNegativeBinomialNLL, subepidemicObjective
from .model import (
    initializeSubepidemicState,
    modifiedLogisticGrowthPatch,
    simulateSubepidemic,
    subepidemicActivationEvents,
    subepidemicIncidence,
    subepidemicOdeOptions,
    subepidemicSimulateRefit,
    unpackSubepidemicParameters,
)
from .noise import AddPoissonError, getMeanVarianceRatio
from .optionsets import (
    GGM,
    GLM,
    GOM,
    GRM,
    LM,
    RICH,
    ForecastOptions,
    Options,
    RtOptions,
    options,
    options_acutelunginjury,
    options_copd,
    options_covid,
    options_covid_wastewater,
    options_diabetes,
    options_ebolaCongo1995,
    options_ebolaUganda2022,
    options_forecast,
    options_forecast_covid,
    options_forecast_ebola_uganda2022,
    options_forecast_wastewater,
    options_measles_SC,
    options_mpx,
    options_Rt,
    options_SARSsingapore2003,
    options_tutorial,
    options_videoYoutube_Faded,
    options_videoYoutube_Uptown_Funk,
    options_wastewater,
    options_yearly,
)
from .plotting import (
    plot_nsubepidemic,
    plot_nsubepidemic_profiles,
    plotFit_subepidemicFramework,
    plotRankings_subepidemicFramework,
    plotReproductionNumber,
)
from .rolling import subepidemicRollingForecast, subepidemicRollingSummary
from .rt import get_Rt
from .runners import Run_Fit_subepidemicFramework
from .scores import (
    computeforecastperformance,
    computeQuantiles,
    computeWIS,
    subepidemicBaselineForecast,
    subepidemicEnsembleWeights,
    subepidemicEvidenceRatios,
    subepidemicWeightedQuantile,
)
from .search import (
    lhsdesign,
    subepidemicOrderLabels,
    subepidemicSearchSpace,
    subepidemicStartPoints,
    subepidemicUnitMap,
)
from .settings import (
    options_fit_settings,
    subepidemicCheckSettings,
    subepidemicFitSeries,
    subepidemicFitSettings,
    subepidemicFitSignature,
    subepidemicSeed,
)
from .utils import (
    Struct,
    datenum,
    datestr,
    datevec,
    example_input_dir,
    get_rng,
    movmean,
    num2ordinal,
    quantile,
    rng,
)

__version__ = "1.0.0"

import types as _types

__all__ = [
    _name
    for _name, _value in list(globals().items())
    if not _name.startswith("_") and not isinstance(_value, _types.ModuleType) and _name != "annotations"
]
