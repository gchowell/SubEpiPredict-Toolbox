"""Options of a run: data set, estimation method, growth model, forecast, Rt.

MATLAB keeps these in ``options.m`` (13 outputs and six global variables),
``options_forecast.m`` and ``options_Rt.m``; an example is run by copying one
of the files in ``examples/`` over ``options.m``. In Python each of them is a
function that returns an object, and the object is passed to the drivers::

    import subepipredict as sp
    sp.Run_Fit_subepidemicFramework(52, '05-11-2020')                    # quick start: sp.options
    sp.Run_Fit_subepidemicFramework(options=sp.options_mpx)              # an example
    sp.Run_Fit_subepidemicFramework(options=sp.options().replace(B=100)) # a changed copy

The rules of ``options.m`` are applied whenever an object is made or changed:
``dist1`` follows ``method1`` for the likelihood methods 1, 3, 4 and 5;
``topmodelsx`` is 1 when ``npatches_fixed`` is 1; with ``onset_fixed = 1`` it
is at most ``npatches_fixed``.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

GGM = 0  # Generalized Growth Model
GLM = 1  # Generalized Logistic Model
GRM = 2  # Generalized Richards Model
LM = 3  # Logistic Model
RICH = 4  # Richards Model
GOM = 5  # Gompertz Model


@dataclass
class Options:
    """The settings of ``options.m``.

    Data
      ``cumulative1``  1 if the input file stores cumulative counts, 0 if incidence
      ``outbreakx``    column of the area or group in the input file (from 1);
                       a list fits several areas in turn
      ``caddate1``     date stamp of the data file, ``'mm-dd-yyyy'``
      ``cadregion``, ``caddisease``, ``datatype``  tags of the file name
      ``DT``           time step: 1 daily, 7 weekly, 365 yearly
      ``datevecfirst1``  date of the first data point ``[yyyy, mm, dd]``
      ``datevecend1``    date of the most recent data file ``[yyyy, mm, dd]``
                         (the file the forecasts are scored against)
    Adjustments
      ``smoothfactor1``       span of the moving average (1 = no smoothing)
      ``calibrationperiod1``  number of most recent points used for calibration
    Estimation
      ``method1``  0 least squares, 1 Poisson likelihood, 3/4/5 negative
                   binomial likelihoods (variance mean+alpha*mean,
                   mean+alpha*mean^2, mean+alpha*mean^d), 6 least absolute deviations
      ``dist1``    error structure of the bootstrap data and the forecast
                   noise: 0 normal, 1 Poisson, 2 negative binomial with
                   VAR = factor*mean, 3/4/5 as the methods, 6 Laplace. Set
                   automatically for methods 1, 3, 4, 5.
      ``numstartpoints``  start points of the multi-start search
      ``B``               bootstrap realizations
    Model
      ``npatches_fixed``  maximum number of sub-epidemics
      ``topmodelsx``      number of top-ranked models kept
      ``flag1``           growth model of every sub-epidemic: 0 GGM, 1 GLM,
                          2 GRM, 3 logistic, 4 Richards, 5 Gompertz
      ``onset_fixed``     1 = all sub-epidemics start together, 0 = asynchronous
    """

    cumulative1: int = 1
    outbreakx: object = 52
    caddate1: str = "05-11-2020"
    cadregion: str = "USA"
    caddisease: str = "coronavirus"
    datatype: str = "deaths"
    DT: int = 1
    datevecfirst1: list = field(default_factory=lambda: [2020, 2, 27])
    datevecend1: list = field(default_factory=lambda: [2022, 5, 9])
    smoothfactor1: int = 7
    calibrationperiod1: int = 90
    method1: int = 3
    dist1: int = 3
    numstartpoints: int = 20
    B: int = 40
    npatches_fixed: int = 2
    topmodelsx: int = 4
    flag1: int = GLM
    onset_fixed: int = 0

    def __post_init__(self):
        self._apply_rules()

    def _apply_rules(self):
        # Automatically set dist1 based on method1 (not for method1 = 0 or 6)
        if self.method1 in (1, 3, 4, 5):
            self.dist1 = self.method1
        # Adjust number of top models if a single subepidemic is used
        if self.npatches_fixed == 1:
            self.topmodelsx = 1
        # Ensure number of top models does not exceed subepidemics when onset is fixed
        try:
            if self.onset_fixed == 1 and self.topmodelsx > self.npatches_fixed:
                self.topmodelsx = self.npatches_fixed
        except TypeError:
            pass  # invalid values are reported by subepidemicCheckSettings
        self.datevecfirst1 = list(self.datevecfirst1)
        self.datevecend1 = list(self.datevecend1)
        for name in ("caddate1", "cadregion", "caddisease", "datatype"):
            value = getattr(self, name)
            if isinstance(value, str):
                setattr(self, name, value.strip())  # strcat drops trailing blanks

    def replace(self, **changes):
        """A copy with some settings changed; the rules of ``options.m`` are applied again."""
        return dataclasses.replace(self, **changes)

    def as_dict(self):
        return dataclasses.asdict(self)

    @property
    def cadtemporal(self):
        return {1: "daily", 7: "weekly", 365: "yearly"}.get(self.DT)


@dataclass
class ForecastOptions:
    """The settings of ``options_forecast.m``.

    ``getperformance``     1 = score the forecast against later data, 0 = skip
    ``deletetempfiles``    1 = delete the temporary ``Forecast-*.mat`` files
    ``forecastingperiod``  forecast horizon (time units ahead)
    ``weight_type1``       ensemble weights: -1 equal, 1 Akaike weights
                           (0 is the same, with a warning), 2 inverse WIS of
                           the calibration period
    """

    getperformance: int = 1
    deletetempfiles: int = 1
    forecastingperiod: int = 4
    weight_type1: int = -1

    def replace(self, **changes):
        return dataclasses.replace(self, **changes)

    def as_dict(self):
        return dataclasses.asdict(self)


@dataclass
class RtOptions:
    """The settings of ``options_Rt.m``.

    ``type_GId1``      generation-interval distribution: 1 gamma, 2 exponential, 3 delta
    ``mean_GI1``       mean of the generation interval (time units of the data)
    ``var_GI1``        variance of the generation interval
    ``minGIcoverage``  Rt is reported only where the curve so far covers at
                       least this share of the generation-interval
                       distribution (0 turns the blanking off)
    """

    type_GId1: int = 1
    mean_GI1: float = 11.7
    var_GI1: float = 3.7
    minGIcoverage: float = 0.99

    def replace(self, **changes):
        return dataclasses.replace(self, **changes)

    def as_dict(self):
        return dataclasses.asdict(self)


def resolve_options(options_=None, default=None):
    """The options object for an ``options`` argument: None (the default
    function), a function returning an object, or the object itself."""
    if options_ is None:
        options_ = default
    if callable(options_):
        options_ = options_()
    return options_


# <============================================================================>
# <============================= options.m (quick start) ======================>
# <============================================================================>


def options():
    """``options.m``: U.S. daily COVID-19 deaths to 05-11-2020 (column 52).

    Negative-binomial likelihood (method1 = 3, dist1 = 3), generalized
    logistic sub-epidemics, up to 2 sub-epidemics, asynchronous onsets,
    7-point moving average for the onset thresholds, calibration period 90,
    20 start points, 4 top models, B = 40. The settings of the README figures.
    """
    return Options()


def options_forecast():
    """``options_forecast.m``: horizon 4, scored, equal ensemble weights."""
    return ForecastOptions()


def options_Rt():
    """``options_Rt.m``: gamma generation interval, mean 11.7, variance 3.7, coverage 0.99."""
    return RtOptions()


# <============================================================================>
# <================================= examples/ ================================>
# <============================================================================>


def options_covid():
    """``examples/options_covid.m``: the quick-start data fitted by least squares."""
    return Options(method1=0, dist1=0)


def options_tutorial():
    """``examples/options_tutorial.m``: layout of the 2024 tutorial paper (least squares, B = 300)."""
    return Options(method1=0, dist1=0, numstartpoints=10, B=300)


def options_SARSsingapore2003():
    """``examples/options_SARSsingapore2003.m``: SARS cases, Singapore 2003 (daily, negative binomial)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="02-25-2003", cadregion="Singapore", caddisease="SARS",
        datatype="cases", DT=1, datevecfirst1=[2003, 2, 25], datevecend1=[2003, 2, 25], smoothfactor1=1,
        calibrationperiod1=70, method1=3, dist1=3, numstartpoints=20, B=300, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_ebolaCongo1995():
    """``examples/options_ebolaCongo1995.m``: Ebola cases, Congo (daily, negative binomial)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="03-08-2023", cadregion="Congo", caddisease="Ebola",
        datatype="cases", DT=1, datevecfirst1=[2022, 8, 1], datevecend1=[2022, 8, 1], smoothfactor1=1,
        calibrationperiod1=129, method1=3, dist1=3, numstartpoints=20, B=300, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_ebolaUganda2022():
    """``examples/options_ebolaUganda2022.m``: Ebola cases, Uganda 2022 (daily, negative binomial)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="11-12-2022", cadregion="Uganda", caddisease="Ebola",
        datatype="cases", DT=1, datevecfirst1=[2022, 8, 8], datevecend1=[2022, 11, 12], smoothfactor1=1,
        calibrationperiod1=97, method1=3, dist1=3, numstartpoints=30, B=300, npatches_fixed=2, topmodelsx=1,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_videoYoutube_Faded():
    """``examples/options_videoYoutube_Faded.m``: YouTube views of *Faded* (daily, least squares)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="01-01-2020", cadregion="World", caddisease="video_Faded",
        datatype="views", DT=1, datevecfirst1=[2020, 1, 1], datevecend1=[2022, 5, 9], smoothfactor1=7,
        calibrationperiod1=86, method1=0, dist1=0, numstartpoints=20, B=300, npatches_fixed=6, topmodelsx=4,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_videoYoutube_Uptown_Funk():
    """``examples/options_videoYoutube_Uptown_Funk.m``: YouTube views of *Uptown Funk* (daily, least squares)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="01-01-2020", cadregion="World", caddisease="video_UptownFunk",
        datatype="views", DT=1, datevecfirst1=[2020, 1, 1], datevecend1=[2022, 5, 9], smoothfactor1=7,
        calibrationperiod1=116, method1=0, dist1=0, numstartpoints=20, B=300, npatches_fixed=6, topmodelsx=4,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_mpx():
    """``examples/options_mpx.m``: U.S. mpox cases to 11-16-2022 (weekly, cumulative file, least squares)."""
    return Options(
        cumulative1=1, outbreakx=1, caddate1="11-16-2022", cadregion="USA", caddisease="monkeypox",
        datatype="cases", DT=7, datevecfirst1=[2022, 5, 10], datevecend1=[2023, 7, 20], smoothfactor1=1,
        calibrationperiod1=10, method1=0, dist1=0, numstartpoints=10, B=300, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_measles_SC():
    """``examples/options_measles_SC.m``: measles cases, South Carolina 2026 (weekly, synchronous onsets)."""
    return Options(
        cumulative1=0, outbreakx=2, caddate1="02-13-2026", cadregion="SC", caddisease="measles",
        datatype="cases", DT=7, datevecfirst1=[2025, 9, 21], datevecend1=[2026, 2, 8], smoothfactor1=1,
        calibrationperiod1=10, method1=0, dist1=0, numstartpoints=20, B=40, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=1,
    )  # fmt: skip


def options_acutelunginjury():
    """``examples/options_acutelunginjury.m``: acute lung injury cases, USA (weekly, least squares)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="03-05-2020", cadregion="USA", caddisease="Acute_lung_injury",
        datatype="cases", DT=7, datevecfirst1=[2020, 3, 31], datevecend1=[2020, 3, 5], smoothfactor1=1,
        calibrationperiod1=46, method1=0, dist1=0, numstartpoints=20, B=300, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_covid_wastewater():
    """``examples/options_covid_wastewater.m``: U.S. COVID-19 wastewater to 12-03-2022 (weekly, column 5)."""
    return Options(
        cumulative1=0, outbreakx=5, caddate1="12-03-2022", cadregion="USA", caddisease="coronavirus",
        datatype="wastewater", DT=7, datevecfirst1=[2022, 1, 1], datevecend1=[2024, 10, 12], smoothfactor1=1,
        calibrationperiod1=10, method1=0, dist1=0, numstartpoints=30, B=300, npatches_fixed=2, topmodelsx=3,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_wastewater():
    """``examples/options_wastewater.m``: U.S. wastewater viral activity level to 01-04-2024 (weekly)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="01-04-2024", cadregion="USA",
        caddisease="Wasterwater_viral_activity_level", datatype="cases", DT=7, datevecfirst1=[2022, 12, 31],
        datevecend1=[2023, 1, 4], smoothfactor1=1, calibrationperiod1=40, method1=0, dist1=0, numstartpoints=20,
        B=300, npatches_fixed=3, topmodelsx=4, flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_copd():
    """``examples/options_copd.m``: COPD deaths, Chatham 2000-2020 (yearly, least squares)."""
    return Options(
        cumulative1=0, outbreakx=1, caddate1="12-31-2020", cadregion="CHATHAM", caddisease="copd",
        datatype="deaths", DT=365, datevecfirst1=[2000, 12, 31], datevecend1=[2020, 12, 31], smoothfactor1=1,
        calibrationperiod1=21, method1=0, dist1=0, numstartpoints=10, B=300, npatches_fixed=2, topmodelsx=4,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_diabetes():
    """``examples/options_diabetes.m``: diabetes deaths, region tag ``usmx``, to 2021 (yearly, column 12)."""
    return Options(
        cumulative1=0, outbreakx=12, caddate1="12-31-2021", cadregion="usmx", caddisease="diabetes",
        datatype="deaths", DT=365, datevecfirst1=[1990, 12, 31], datevecend1=[2021, 12, 31], smoothfactor1=7,
        calibrationperiod1=32, method1=0, dist1=0, numstartpoints=20, B=300, npatches_fixed=3, topmodelsx=4,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_yearly():
    """``examples/options_yearly.m``: the mpox series relabeled as yearly (placeholder dates)."""
    return Options(
        cumulative1=1, outbreakx=1, caddate1="11-16-2049", cadregion="USA", caddisease="monkeypox",
        datatype="cases", DT=365, datevecfirst1=[2022, 5, 10], datevecend1=[2079, 7, 20], smoothfactor1=1,
        calibrationperiod1=10, method1=0, dist1=0, numstartpoints=10, B=300, npatches_fixed=2, topmodelsx=2,
        flag1=GLM, onset_fixed=0,
    )  # fmt: skip


def options_forecast_covid():
    """``examples/options_forecast_covid.m``: horizon 30, scored, Akaike weights."""
    return ForecastOptions(getperformance=1, deletetempfiles=1, forecastingperiod=30, weight_type1=1)


def options_forecast_ebola_uganda2022():
    """``examples/options_forecast_ebola_uganda2022.m``: horizon 30, not scored, equal weights."""
    return ForecastOptions(getperformance=0, deletetempfiles=1, forecastingperiod=30, weight_type1=-1)


def options_forecast_wastewater():
    """``examples/options_forecast_wastewater.m``: horizon 4, not scored, equal weights."""
    return ForecastOptions(getperformance=0, deletetempfiles=1, forecastingperiod=4, weight_type1=-1)
