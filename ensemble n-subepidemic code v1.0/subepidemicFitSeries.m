function [ydata,info] = subepidemicFitSeries(data,smoothfactor1,method1,settings)
%SUBEPIDEMICFITSERIES The series the objective is fitted to.
%   [ydata,info] = subepidemicFitSeries(data,smoothfactor1,method1,settings)
%
%   The screen and the bootstrap both call this function, so they cannot
%   disagree about what was fitted.
%
%   Least squares (method 0) and absolute deviations (method 6):
%       ydata = moving average of DATA with span smoothfactor1, as before.
%
%   Count likelihoods (methods 1, 3, 4, 5), default:
%       ydata = DATA, the raw counts, whatever smoothfactor1 is.
%       smoothfactor1 is still used by the screen to build the grid of onset
%       thresholds from the smoothed cumulative curve; it no longer touches
%       the likelihood.
%
%   Why. A Poisson or negative-binomial likelihood describes the scatter of
%   the counts about the mean curve. A moving average removes most of that
%   scatter before the likelihood sees it, so the fitted dispersion is too
%   small and prediction intervals are too narrow for the observed counts.
%   In the bundled COVID-19 run (method 3, smoothfactor1 = 7) the fit to the
%   7-day average gave alpha = 1.93; its 95% interval covered 56% of the
%   calibration counts and 20% of the next 30 days.
%
%   settings.countLikelihoodOnRawData = false restores the previous
%   behaviour (count likelihoods fitted to the moving average); the screen
%   then warns.
%
%   Count likelihoods need nonnegative data. A negative value in the raw
%   series (for example a negative increment of a corrected cumulative
%   series) is an error here, with the choices spelled out, because every
%   negative-binomial fit would otherwise fail without a clear reason.
%   Non-integer raw values are accepted with a warning: the toolbox's
%   likelihoods keep their existing convention for fractional data.
%
%   info  struct: usedRawData (logical), smoothingSpan (the span applied to
%         ydata, 1 when none), smoothedCountLikelihood (logical: a count
%         likelihood is being fitted to smoothed data), nNonInteger.
%
%   See also subepidemicFitSettings, fittingModifiedLogisticFunctionPatchABC

if nargin<4 || isempty(settings)
    settings = subepidemicFitSettings();
end
if ~(isnumeric(data) && isreal(data) && isvector(data) && ~isempty(data) && all(isfinite(data(:))))
    error('SubEpiPredict:InvalidData','data must be a nonempty finite real vector.');
end
if ~(isnumeric(method1) && isscalar(method1) && any(method1==0:6))
    error('SubEpiPredict:InvalidMethod','method1 must be an integer from 0 to 6.');
end
if isempty(smoothfactor1)
    smoothfactor1 = 1;
end
if ~(isnumeric(smoothfactor1) && isreal(smoothfactor1) && isscalar(smoothfactor1) && ...
        isfinite(smoothfactor1) && smoothfactor1>=1)
    error('SubEpiPredict:InvalidSmoothing','smoothfactor1 must be a number of at least 1.');
end
data = double(data(:));

isCount = any(method1==[1 3 4 5]);
useRaw = isCount && settings.countLikelihoodOnRawData;

info = struct('usedRawData',false,'smoothingSpan',1, ...
    'smoothedCountLikelihood',false,'nNonInteger',0);

if useRaw || smoothfactor1<=1
    ydata = data;
    info.usedRawData = true;
else
    ydata = smoothdata(data,'movmean',smoothfactor1);
    info.smoothingSpan = smoothfactor1;
    info.smoothedCountLikelihood = isCount;
end

if isCount && info.usedRawData
    nNegative = sum(data<0);
    if nNegative>0
        error('SubEpiPredict:NegativeCounts', ...
            ['The calibration series has %d negative value(s), and method1 = %d fits a count likelihood ', ...
             'to the raw series. Correct the data, or use least squares (method1 = 0), or set ', ...
             'countLikelihoodOnRawData = false in options_fit_settings.m to fit the moving average as before.'], ...
            nNegative,method1);
    end
    info.nNonInteger = sum(data~=round(data));
    if info.nNonInteger>0
        subepidemicWarning('SubEpiPredict:NonIntegerCounts', ...
            ['The calibration series has %d non-integer value(s) and method1 = %d fits a count likelihood. ', ...
             'The likelihood keeps its convention for fractional data; integer counts are what it describes.'], ...
            info.nNonInteger,method1);
    end
end
end
