function U = options_fit_settings_example()
%OPTIONS_FIT_SETTINGS_EXAMPLE Template for optional fit settings.
%   Copy this file to options_fit_settings.m (and rename the function) to
%   override defaults of subepidemicFitSettings.m. Return only the fields
%   you want to change; every other setting keeps its default. While no
%   options_fit_settings.m exists, the defaults are used.

U = struct();

% ---- bounds (empty or absent = default bound, see subepidemicSearchSpace.m) --
% The lower bound of r is the legacy one, mean(first two observations)/200.
% In tests with calibration windows that start mid-epidemic the fits stayed
% above it (they use a large r with a small p), so it was kept. A fixed floor
% instead:
% U.rLB = 1e-3;
%
% The upper bound of r is the larger of the legacy bound, 50*max(first two
% observations), and 10 times the largest count. A fixed value instead:
% U.rUB = 500;
% The legacy bound alone, as before (50 for a series that starts with 1, 0):
% U.rUpperBoundFactor = 0;
%
% Carrying capacity between 20 and 100 times the cumulative count, say 5e6:
% U.KUB = 5e6;
%
% Fix p of every sub-epidemic at 1 (lower bound = upper bound):
% U.pLB = 1;  U.pUB = 1;
%
% One value per sub-epidemic:
% U.KUB = [2e6 5e6];

% ---- speed ------------------------------------------------------------------
% Local solves on a parallel pool (Parallel Computing Toolbox):
% U.useParallel = true;
%
% Looser fmincon tolerances (faster, slightly less precise optima):
% U.optimOptions = struct('StepTolerance',1e-8,'FunctionTolerance',1e-8);

% ---- behaviour --------------------------------------------------------------
% Fit Poisson and negative-binomial likelihoods to the moving average of span
% smoothfactor1, as the toolbox did before (default: the raw counts):
% U.countLikelihoodOnRawData = false;
%
% Fit every threshold from the common start set only:
% U.continuationStart = false;
%
% More exploration in each bootstrap refit:
% U.bootstrapStarts = 3;

% ---- repeatable runs ----------------------------------------------------------
% The random number generator is reset with this seed at the start of the fit
% (once per area) and of the forecast, fit-plot and Rt scripts, so the same
% call on the same files gives the same numbers again. The default is 1 since
% revision l (it was empty before). Empty leaves the generator alone, so every
% run draws new bootstrap data and forecast noise; the fit then still saves
% the generator state with its results (see subepidemicSeed.m):
% U.rngSeed = [];

% ---- refits that cannot be simulated in a forecast ------------------------------
% The forecast and Rt scripts extend every bootstrap refit past the
% calibration window. A refit whose curve cannot be integrated that far is
% left out with a warning, as long as the share left out for a ranked model
% stays at or below this value (default 0.05); above it the script stops.
% 0 stops at the first such refit; 1 stops only when none is left:
% U.maxSkippedRefitFraction = 0;

% ---- baseline of the forecast scores --------------------------------------------
% With getperformance = 1 the forecast script also scores a flat baseline: the
% mean of the last observations, repeated over the forecast period
% (subepidemicBaselineForecast.m). Empty (the default) takes the last 7
% observations of daily data (DT = 1) and the last one otherwise. Another
% window:
% U.baselineWindow = 14;
% ---- identifiability report (PII) ------------------------------------------------
% The fit and plotFit_subepidemicFramework report, for every parameter, the
% parameter identifiability index PII = log10(UB/LB) of its 95% bootstrap
% interval and its companions (subepidemicIdentifiability.m). The headline
% measure is the PII where the lower search bound is positive (r, K, alpha, d)
% and the share of the search range otherwise (p, and a of the GRM and Richards
% kernels). To choose it yourself give one entry per parameter type,
% [r p a K alpha d]: 1 = log scale, 0 = range fraction, NaN = automatic.
% Range fraction for r, automatic for the rest:
% U.PIIscale = [0 NaN NaN NaN NaN NaN];

% ---- weakly supported top models --------------------------------------------------
% The fit reports the evidence ratio of each top-ranked model against ranked
% model 1, exp(dAICc/2), and warns (SubEpiPredict:WeakTopModel) when it is at
% least this value; default 100, an AICc difference of 9.2. Report only: the
% top models are kept. Inf turns the warning off:
% U.evidenceRatioLimit = 1000;

% ---- start points of the first threshold --------------------------------------------
% A fit of the screen without a solution of a previous onset threshold to
% start from (the first threshold of each sweep) gets this many times
% numstartpoints start points; default 3. 1 gives no extra start points:
% U.firstThresholdStartFactor = 1;

% ---- parameters that are not identified ------------------------------------------
% A parameter whose 95% bootstrap interval covers at least this share of what
% the search bounds allow is reported as not identified (warning
% SubEpiPredict:NotIdentified); default 0.5. Inf turns the warning off:
% U.notIdentifiedRatio = 0.8;
end
