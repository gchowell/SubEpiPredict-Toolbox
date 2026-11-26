# Input data

Each file is a header-free numeric matrix: one row per regularly spaced observation, one column per area or group. The fit and forecast drivers read files from this folder by name, so the name must match the settings in `options.m`.

## Naming convention

```
[cumulative-]<frequency>-<caddisease>-<datatype>-<cadregion>-<mm-dd-yyyy>.txt
```

- `cumulative-` is present when `cumulative1 = 1` (running totals) and absent for interval counts.
- `<frequency>` is `daily`, `weekly` or `yearly` (`DT = 1`, `7` or `365`).
- The date is always written with a four-digit year. The fitting file uses `caddate1`; the evaluation file uses `datevecend1`.
- Names are case-sensitive on Linux, so `caddisease`, `datatype` and `cadregion` must match the filename exactly.

## Bundled datasets

| Files | Used by |
|---|---|
| `cumulative-daily-coronavirus-deaths-USA-<date>.txt` (snapshots taken Apr–May 2020) and `...-05-09-2022.txt` (later series for scoring) | `options.m` (quick start), `options_covid.m`, `options_subepi.m`, `options_tutorial.m`, rolling-origin example |
| `cumulative-daily-coronavirus-cases-USA-<date>.txt` (Apr–May 2020 snapshots) | Available for case forecasts |
| `cumulative-weekly-monkeypox-cases-USA-11-16-2022.txt`, `...-07-20-2023.txt` | `options_mpx.m`, `options_mpx1.m` |
| `cumulative-yearly-monkeypox-cases-USA-11-16-2049.txt`, `...-07-20-2079.txt` | `options_yearly.m` (the weekly series relabeled to demonstrate `DT = 365`; the dates are placeholders) |
| `daily-SARS-cases-Singapore-02-25-2003.txt` | `options_SARSsingapore2003.m` |
| `daily-Ebola-cases-Congo-03-08-2023.txt` | `options_ebolaCongo1995.m` |
| `daily-Ebola-cases-Uganda-11-12-2022.txt` | `options_ebolaUganda2022.m` |
| `daily-video_UptownFunk-views-World-01-01-2020.txt`, `daily-video_Faded-views-World-01-01-2020.txt` | `options_videoYoutube*.m` |
| `weekly-Acute_lung_injury-cases-USA-03-05-2020.txt` | `options_acutelunginjury.m` |
| `weekly-coronavirus-wastewater-USA-<date>.txt` | `options_faharudeen.m` |
| `weekly-Wasterwater_viral_activity_level-cases-USA-01-04-2024.txt` | `options_wastewater.m` |
| `weekly-measles-cases-SC-02-08-2026.txt`, `...-02-13-2026.txt` | `options_measles_SC.m` |
| `yearly-copd-deaths-CHATHAM-12-31-2020.txt` | `options_copd.m` |
| `yearly-diabetes-deaths-usmx-12-31-2021.txt` | `options_diabetes.m` |
| `yearly-HIV-cases-World_Regions-01-01-2019.txt` | Available as an additional example |
| `cumulative-daily-coronavirus-deaths-USA-04-27-20.txt` | Original real-time snapshot of 04-27-2020 (see below) |
| `Cumulative_Most_Recent_Timeseries_US-CDC.txt` | Not read by the toolbox (two columns: day index and cumulative count) |

`cumulative-daily-coronavirus-deaths-USA-04-27-2020.txt` matches the later 2022 series through April 27 in every column. `...-04-27-20.txt` is the data as reported at the time: column 52 (USA) is identical, but 12 state columns differ because of later reporting revisions.
