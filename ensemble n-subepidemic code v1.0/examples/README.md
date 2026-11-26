# Example configurations

The drivers always read `options.m` and `options_forecast.m` in the code folder. The files here are ready-made settings for the bundled datasets. They are kept in this subfolder so they cannot be picked up by mistake.

**To run an example,** copy the contents of one `options_<name>.m` file over `options.m` in the code folder (keep the function name `options`). If the example has a matching `options_forecast_<name>.m`, copy that over `options_forecast.m` as well. Then call `Run_Fit_subepidemicFramework` and `plotForecast_subepidemicFramework` as described in the main README. Keep a copy of the original `options.m` to return to the quick start.

## Fit settings (`options.m`)

| File | Data in `input/` | Frequency | Calibration rows | Fitting method |
|---|---|---|---:|---|
| `options_covid.m` | U.S. COVID-19 deaths to 05-11-2020 | daily | 90 | least squares (quick start uses negative binomial) |
| `options_tutorial.m` | U.S. COVID-19 deaths to 05-11-2020 | daily | 90 | least squares, layout of the 2024 tutorial paper |
| `options_SARSsingapore2003.m` | SARS cases, Singapore 2003 | daily | 70 | negative binomial |
| `options_ebolaCongo1995.m` | Ebola cases, Congo (1995 outbreak) | daily | 129 | negative binomial |
| `options_ebolaUganda2022.m` | Ebola cases, Uganda 2022 | daily | 97 | negative binomial |
| `options_videoYoutube_Faded.m` | YouTube views, *Faded* | daily | 86 | least squares |
| `options_videoYoutube_Uptown_Funk.m` | YouTube views, *Uptown Funk* | daily | 116 | least squares |
| `options_mpx.m` | U.S. mpox cases to 11-16-2022 | weekly | 10 | least squares |
| `options_measles_SC.m` | Measles cases, South Carolina 2026 | weekly | 10 | least squares |
| `options_acutelunginjury.m` | Acute lung injury cases, USA | weekly | 46 | least squares |
| `options_covid_wastewater.m` | U.S. COVID-19 wastewater to 12-03-2022 | weekly | 10 | least squares |
| `options_wastewater.m` | U.S. wastewater viral activity level to 01-04-2024 | weekly | 40 | least squares |
| `options_copd.m` | COPD deaths, Chatham 2000–2020 | yearly | 21 | least squares |
| `options_diabetes.m` | Diabetes deaths, region tag `usmx`, to 2021 | yearly | 32 | least squares |
| `options_yearly.m` | Mpox series relabeled as yearly (placeholder dates) | yearly | 10 | least squares |

## Forecast settings (`options_forecast.m`)

| File | Horizon | Scores the forecast (`getperformance`) | Ensemble weights |
|---|---:|---|---|
| `options_forecast_covid.m` | 30 | yes | Akaike (`1`) |
| `options_forecast_ebola_uganda2022.m` | 30 | no | equal (`-1`) |
| `options_forecast_wastewater.m` | 4 | no | equal (`-1`) |

Scoring a forecast needs a later data file named by `datevecend1`. The Congo Ebola, YouTube and wastewater-viral-activity examples have no such file in `input/`, so run them with `getperformance = 0`.
