# Legacy code

These functions come from earlier versions of the toolbox. None of the current drivers calls them, and they are not on the MATLAB path when you run the toolbox from the code folder. They are kept for reference and for older scripts.

| File | What it was |
|---|---|
| `getbounds.m` | Parameter bounds, now built by `subepidemicSearchSpace.m` |
| `getensemblesubepidemics_Deaths.m` | Earlier ensemble forecast routine; reads output filenames that the current drivers no longer write |
| `plotModifiedLogisticGrowthPatchMethodLogLik.m`, `plotModifiedLogisticGrowthPatchMethodsFullLogLik.m` | Earlier objective functions that read global variables; replaced by `subepidemicObjective.m` |
| `sesweights.m` | Exponential smoothing weights, not used by the current ensemble weighting |

To call one from an old script, run `addpath('legacy')` from the code folder first.
