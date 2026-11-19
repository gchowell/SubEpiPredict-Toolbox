function S = subepidemicFitSettings()
%SUBEPIDEMICFITSETTINGS Optimizer and search-space settings of the fitting layer.
%   S = subepidemicFitSettings() returns the defaults below. If a function
%   options_fit_settings.m is on the path, it is called with no inputs and
%   must return a struct; each of its fields replaces the default of the same
%   name. Unknown field names are an error, so a misspelled setting cannot be
%   ignored silently. Copy options_fit_settings_example.m to
%   options_fit_settings.m to start. Existing options*.m files need no change.
%
%   SEARCH SPACE (subepidemicSearchSpace.m)
%     logScaleR, logScaleK, logScaleAlpha   true: r, K and the negative-
%         binomial alpha are searched on a log scale. The bounds themselves do
%         not change; only the coordinates the optimizer and the Latin
%         hypercube work in. With K in [20, 1e10] a linear scale puts 20 of 21
%         start points above 10 times the total count.
%     rLB, rUB, pLB, pUB, aLB, aUB, KLB, KUB   [] = the legacy bound of
%         getbounds.m. A scalar applies to every sub-epidemic; a vector gives
%         one value per sub-epidemic (the first n entries are used). Setting a
%         lower bound equal to its upper bound fixes the parameter.
%     rUpperBoundFactor   10: the default upper bound of r is the larger of
%         the legacy bound, 50*max(first two observations), and 10 times the
%         largest count of the calibration series. The legacy bound depends
%         on how the series happens to start (50 when it starts with 1, 0)
%         and was reached by fits and bootstrap refits. 0 gives the legacy
%         bound alone. Has no effect where rUB is set.
%
%   FITTED SERIES (subepidemicFitSeries.m)
%     countLikelihoodOnRawData   true: Poisson and negative-binomial fits
%         (methods 1, 3, 4, 5) use the raw counts whatever smoothfactor1 is;
%         smoothfactor1 then only shapes the grid of onset thresholds.
%         false: they are fitted to the moving average, as before (the screen
%         warns). Least squares and absolute deviations are not affected.
%
%   OPTIMIZER (fitSubepidemicModel.m)
%     optimOptions      struct of fmincon option name/value pairs applied on
%                       top of the defaults (sqp, StepTolerance 1e-10,
%                       FunctionTolerance 1e-10, OptimalityTolerance 1e-8,
%                       forward differences with step 1e-6 in unit coordinates).
%     fitMaxAttempts    3: MultiStart passes with fresh start points before a
%                       candidate is reported as failed.
%     useParallel       false. true runs the local solves of each MultiStart
%                       call on a parallel pool (Parallel Computing Toolbox).
%     polish            true: Nelder-Mead pass after MultiStart, run only when
%                       the stationarity check fails.
%     polishGradientTol 10: the polish runs when the largest gradient
%                       component exceeds this value. Units: log-likelihood
%                       units per unit coordinate (a unit coordinate is the
%                       full range of a parameter).
%     polishBudget      100: Nelder-Mead evaluations per free parameter.
%     agreementTol      0.01: two local runs agree when their objectives differ
%                       by less than this many log-likelihood units.
%
%   SCREEN AND BOOTSTRAP
%     continuationStart true: in the threshold screen, the accepted solution
%                       of the previous threshold (same number of
%                       sub-epidemics) is added to the start points of the next.
%     bootstrapStarts   1: Latin-hypercube start points per bootstrap refit, in
%                       addition to the fitted vector and one jittered copy of it.
%
%   REPEATABLE RUNS (subepidemicSeed.m)
%     rngSeed           1: the random number generator is reset with this
%                       seed at the start of the fit (once per area) and of
%                       the forecast, fit-plot and Rt scripts, so the same call
%                       on the same files gives the same numbers. [] leaves
%                       the generator alone: every run draws new bootstrap
%                       data and forecast noise. Until revision l the default
%                       was [].
%
%   SCREEN (fittingModifiedLogisticFunctionPatchABC.m)
%     firstThresholdStartFactor 3: a fit of the screen that has no solution of
%                       a previous onset threshold to start from (the first
%                       threshold of each sweep) gets this many times
%                       numstartpoints start points. 1: no extra start points,
%                       as until revision l.
%
%   IDENTIFIABILITY WARNING (fittingModifiedLogisticFunctionPatchMultiple.m)
%     notIdentifiedRatio 0.5: a parameter whose 95% bootstrap interval covers
%                       at least this share of what the search bounds allow
%                       (PII_ratio of the PII report) is reported as not
%                       identified (warning SubEpiPredict:NotIdentified).
%                       Inf turns the warning off.
%
%   RANKING (subepidemicEvidenceRatios.m)
%     evidenceRatioLimit 100: a top-ranked model whose evidence ratio against
%                       ranked model 1, exp(dAICc/2), is at least this value
%                       is reported as weakly supported (warning
%                       SubEpiPredict:WeakTopModel). 100 is an AICc difference
%                       of 9.2. Report only: the top models are kept. Inf
%                       turns the warning off.
%
%   IDENTIFIABILITY REPORT (subepidemicIdentifiability.m)
%     PIIscale          []: the headline measure of each parameter is chosen
%                       automatically, the PII = log10(UB/LB) of the 95%
%                       interval when the lower search bound is positive (r, K,
%                       alpha, d) and the range fraction otherwise (p, and a
%                       of the GRM and Richards kernels). A vector of 6
%                       entries, one per parameter type [r p a K alpha d],
%                       overrides it: 1 = log scale, 0 = range fraction,
%                       NaN = automatic. It changes PII_scale and PII_ratio
%                       of the report only, never the fit.
%
%   See also fitSubepidemicModel, subepidemicSearchSpace, subepidemicOdeOptions

S = struct( ...
    'logScaleR',true,'logScaleK',true,'logScaleAlpha',true, ...
    'rLB',[],'rUB',[],'pLB',[],'pUB',[],'aLB',[],'aUB',[],'KLB',[],'KUB',[], ...
    'rUpperBoundFactor',10, ...
    'countLikelihoodOnRawData',true, ...
    'optimOptions',struct(), ...
    'fitMaxAttempts',3, ...
    'useParallel',false, ...
    'polish',true,'polishGradientTol',10,'polishBudget',100, ...
    'agreementTol',1e-2, ...
    'continuationStart',true, ...
    'bootstrapStarts',1, ...
    'rngSeed',1, ...
    'maxSkippedRefitFraction',0.05, ...
    'baselineWindow',[], ...
    'PIIscale',[], ...
    'evidenceRatioLimit',100, ...
    'firstThresholdStartFactor',3, ...
    'notIdentifiedRatio',0.5);

if exist('options_fit_settings','file')==2
    U = options_fit_settings();
    if ~(isstruct(U) && isscalar(U))
        error('SubEpiPredict:FitSettings', ...
            'options_fit_settings must return one struct of settings.');
    end
    names = fieldnames(U);
    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('SubEpiPredict:FitSettings', ...
                'Unknown fit setting ''%s''. See subepidemicFitSettings.m for the list.',names{k});
        end
        S.(names{k}) = U.(names{k});
    end
end

% ---- validation -------------------------------------------------------------
flags = {'logScaleR','logScaleK','logScaleAlpha','countLikelihoodOnRawData','useParallel','polish','continuationStart'};
for k = 1:numel(flags)
    v = S.(flags{k});
    if ~((islogical(v) || isnumeric(v)) && isscalar(v) && (v==0 || v==1))
        error('SubEpiPredict:FitSettings','Setting %s must be true or false.',flags{k});
    end
    S.(flags{k}) = logical(v);
end
positive = {'polishGradientTol','polishBudget','agreementTol'};
for k = 1:numel(positive)
    v = S.(positive{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>0)
        error('SubEpiPredict:FitSettings','Setting %s must be a positive number.',positive{k});
    end
end
v = S.rUpperBoundFactor;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0)
    error('SubEpiPredict:FitSettings','Setting rUpperBoundFactor must be a number of at least 0.');
end
counts = {'fitMaxAttempts','bootstrapStarts'};
minimum = [1 0];
for k = 1:numel(counts)
    v = S.(counts{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v==fix(v) && v>=minimum(k))
        error('SubEpiPredict:FitSettings','Setting %s must be an integer of at least %d.',counts{k},minimum(k));
    end
end
bounds = {'rLB','rUB','pLB','pUB','aLB','aUB','KLB','KUB'};
for k = 1:numel(bounds)
    v = S.(bounds{k});
    if ~(isempty(v) || (isnumeric(v) && isreal(v) && isvector(v) && all(isfinite(v(:)))))
        error('SubEpiPredict:FitSettings', ...
            'Setting %s must be empty, a finite scalar, or a finite vector.',bounds{k});
    end
end
v = S.rngSeed;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0 && v==fix(v) && v<2^32))
    error('SubEpiPredict:FitSettings','Setting rngSeed must be empty or a nonnegative integer below 2^32.');
end
v = S.baselineWindow;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=1 && v==fix(v)))
    error('SubEpiPredict:FitSettings','Setting baselineWindow must be empty or a positive integer.');
end
v = S.maxSkippedRefitFraction;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && v>=0 && v<=1)
    error('SubEpiPredict:FitSettings','Setting maxSkippedRefitFraction must be a number from 0 to 1.');
end
v = S.firstThresholdStartFactor;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=1)
    error('SubEpiPredict:FitSettings','Setting firstThresholdStartFactor must be a number of at least 1 (1: no extra start points).');
end
v = S.notIdentifiedRatio;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && ~isnan(v) && v>0)
    error('SubEpiPredict:FitSettings','Setting notIdentifiedRatio must be a positive number (Inf turns the warning off).');
end
v = S.evidenceRatioLimit;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && ~isnan(v) && v>1)
    error('SubEpiPredict:FitSettings','Setting evidenceRatioLimit must be a number above 1 (Inf turns the warning off).');
end
v = S.PIIscale;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && numel(v)==6 && all(isnan(v(:)) | v(:)==0 | v(:)==1)))
    error('SubEpiPredict:FitSettings', ...
        'Setting PIIscale must be empty or 6 entries [r p a K alpha d], each 1 (log scale), 0 (range fraction) or NaN (automatic).');
end
if ~(isstruct(S.optimOptions) && isscalar(S.optimOptions))
    error('SubEpiPredict:FitSettings','Setting optimOptions must be a struct of fmincon options.');
end
end
