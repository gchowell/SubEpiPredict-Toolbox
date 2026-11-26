function B = subepidemicBaselineForecast(y,horizon,window,nSamples)
%SUBEPIDEMICBASELINEFORECAST Flat baseline forecast with empirical uncertainty.
%   B = subepidemicBaselineForecast(y,horizon,window) forecasts HORIZON steps
%   beyond the series Y (the calibration data, one value per time unit) without
%   any model. It is the reference a model forecast should beat.
%
%   Point forecast. The mean of the last WINDOW observations, for every
%   horizon. WINDOW = 7 smooths the weekly reporting cycle of daily data;
%   WINDOW = 1 is the last observation (weekly or yearly data).
%
%   Uncertainty. For horizon h, the same rule is applied at every earlier
%   time t of the series, and its error is recorded:
%       e = y(t+h) - mean(y(t-window+1 : t)),   t = window, ..., n-h.
%   The forecast distribution at horizon h is the point forecast plus these
%   errors and their negatives. Using both signs keeps the median at the point
%   forecast, as in the baseline of the COVID-19 Forecast Hub; without it the
%   rise of the epidemic so far would be projected forward. Values below 0
%   are set to 0. If fewer than 5 errors exist for a horizon, the errors of
%   the longest shorter horizon that has 5 are used; if none has, the
%   uncertainty is not defined (NaN) and only point scores can be computed.
%
%   B.level      the point forecast (one number)
%   B.point      horizon-by-1, all equal to B.level
%   B.samples    horizon-by-nSamples (default 200): evenly spaced quantiles of
%                the forecast distribution of each horizon, so that quantiles
%                of a row reproduce that distribution. No random numbers.
%   B.nErrors    horizon-by-1, number of past errors behind each horizon
%   B.borrowed   horizon-by-1, true where a shorter horizon's errors were used
%   B.window     the window used
%
%   The forecast script scores B with the functions that score the models
%   (computeWIS.m, computeforecastperformance.m).

if nargin<4 || isempty(nSamples)
    nSamples = 200;
end
if ~(isnumeric(y) && isreal(y) && isvector(y) && ~isempty(y) && all(isfinite(y(:))))
    error('SubEpiPredict:Baseline','The series must be a non-empty vector of finite numbers.');
end
if ~(isnumeric(horizon) && isscalar(horizon) && isfinite(horizon) && horizon>=1 && horizon==fix(horizon))
    error('SubEpiPredict:Baseline','The horizon must be a positive integer.');
end
if ~(isnumeric(window) && isscalar(window) && isfinite(window) && window>=1 && window==fix(window))
    error('SubEpiPredict:Baseline','The window must be a positive integer.');
end
y = double(y(:));
n = numel(y);
window = min(window,n);

B = struct('level',mean(y(n-window+1:n)),'point',[],'samples',NaN(horizon,nSamples), ...
    'nErrors',zeros(horizon,1),'borrowed',false(horizon,1),'window',window);
B.point = repmat(B.level,horizon,1);

% rolling mean of the last `window` values, available from t = window on
c = cumsum([0;y]);
level = (c(window+1:n+1)-c(1:n-window+1))/window;   % level(k) is the mean ending at t = window+k-1

minErrors = 5;
probabilities = ((1:nSamples)-0.5)/nSamples;
lastGood = [];
for h = 1:horizon
    nOrigins = n-h-window+1;                         % origins t = window, ..., n-h
    if nOrigins>=1
        e = y(window+h:n)-level(1:nOrigins);
    else
        e = [];
    end
    B.nErrors(h) = numel(e);
    if numel(e)>=minErrors
        lastGood = e;
    elseif ~isempty(lastGood)
        e = lastGood;
        B.borrowed(h) = true;
        B.nErrors(h) = numel(e);
    else
        continue                                     % uncertainty not defined at this horizon
    end
    both = sort([e;-e]);
    B.samples(h,:) = max(B.level+localQuantile(both,probabilities),0);
end
end

% -----------------------------------------------------------------------------
function q = localQuantile(x,p)
% Quantiles of the sorted vector x by the rule of MATLAB's quantile: linear
% interpolation with the k-th value at probability (k-0.5)/n.
n = numel(x);
position = min(max(p*n+0.5,1),n);
lower = floor(position);
upper = min(lower+1,n);
weight = position-lower;
q = x(lower)'.*(1-weight)+x(upper)'.*weight;
end
