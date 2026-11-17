function [Rs,ts,coverage]=get_Rt(timevect,incidence1,type_GId1,mean1,var1,minCoverage)
%GET_RT Effective reproduction number from an incidence curve (renewal equation).
%   [Rs,ts] = get_Rt(timevect,incidence1,type_GId1,mean1,var1)
%   [Rs,ts,coverage] = get_Rt(timevect,incidence1,type_GId1,mean1,var1,minCoverage)
%
%   timevect, incidence1   times and model incidence (same length)
%   type_GId1              generation-interval distribution:
%                          1 = gamma, 2 = exponential, 3 = delta (fixed interval)
%   mean1, var1            mean and VARIANCE of the generation interval, in the
%                          time units of timevect (var1 is used by the gamma only)
%   minCoverage            optional, default 0.99. Rt is reported only where the
%                          curve so far covers at least this share of the
%                          generation-interval distribution; before that it is
%                          NaN. Use 0 for the earlier behaviour (no blanking).
%
%   Rs(k) is the estimate at time ts(k+1); ts is the grid of step 0.01 from
%   timevect(1) to timevect(end). coverage(k) is the share of the
%   generation-interval distribution that lies within the curve at ts(k+1).
%
%   Why early values are blanked. The estimate is
%       R(t) = I(t) / sum over past lags s of I(t-s) w(s),
%   and the sum can only run over the part of the curve that exists. At a
%   time t after the start of the curve, only lags up to t are available, so
%   the denominator misses the share 1 - coverage of the infectors and R(t)
%   is too large by a factor of about 1/coverage. For a constant incidence
%   and the default gamma interval (mean 11.7, variance 3.7) the unblanked
%   estimate is in the hundreds of thousands at day 5 and 5.2 at day 10,
%   where the true value is 1. From coverage 0.99 on, that error is at most
%   about 1 percent.
%
%   Changes from the earlier get_Rt.m: the blanking above; the delta interval
%   is placed at lag mean1 exactly (it was one grid step late, and a mean
%   that is not a multiple of 0.01 in floating point stopped with an indexing
%   error); the inputs are checked.

if nargin<6 || isempty(minCoverage)
    minCoverage=0.99;
end
if nargin<5
    var1=[];
end
if ~(isnumeric(minCoverage) && isreal(minCoverage) && isscalar(minCoverage) && minCoverage>=0 && minCoverage<=1)
    error('SubEpiPredict:InvalidRtSettings','The minimum generation-interval coverage must be one number from 0 to 1.');
end
if ~(isnumeric(type_GId1) && isscalar(type_GId1) && any(type_GId1==[1 2 3]))
    error('SubEpiPredict:InvalidRtSettings','type_GId1 must be 1 (gamma), 2 (exponential) or 3 (delta).');
end
if ~(isnumeric(mean1) && isreal(mean1) && isscalar(mean1) && isfinite(mean1) && mean1>0)
    error('SubEpiPredict:InvalidRtSettings','The mean of the generation interval must be one finite positive number.');
end
if type_GId1==1 && ~(isnumeric(var1) && isreal(var1) && isscalar(var1) && isfinite(var1) && var1>0)
    error('SubEpiPredict:InvalidRtSettings','The gamma generation interval needs a finite positive variance.');
end
if numel(timevect)~=numel(incidence1)
    error('SubEpiPredict:InvalidRtSettings','timevect and the incidence curve must have the same length.');
end

if numel(timevect)<2
    Rs = []; ts = timevect; coverage = [];   % one point: nothing to estimate
    return
end

dt1=0.01;
timevect2=timevect(1):dt1:timevect(end);

incidence1=interp1(timevect,incidence1,timevect2);

timevect3=timevect2-min(timevect2);

ts=timevect2;

n = numel(incidence1);
if n < 2
    Rs = []; coverage = [];
    return
end

% Generation-interval weights: prob1(j) is the probability of a lag in
% ((j-1)*dt1, j*dt1], used at lag j*dt1.

if type_GId1==1
    % gamma
    b=var1/mean1;
    a=mean1/b;
    prob1=diff(gamcdf(timevect3,a,b))';
elseif type_GId1==2
    % exponential
    prob1=diff(expcdf(timevect3,mean1))';
else
    % delta: all weight at lag mean1 (to the nearest grid step, at least one)
    lag1=max(1,round(mean1/dt1));
    prob1=zeros(max(n-1,lag1),1);
    prob1(lag1)=1;
end

% Renewal denominator by full convolution: at time index i it is
% c(i-1) = sum_{j=1}^{i-1} incidence1(i-j) * prob1(j).

I = incidence1(:);
w = prob1(:);
c = conv(I, w, 'full');
denominator = c(1:n-1);

Rs = NaN(n-1, 1);
currentIncidence = I(2:n);
valid = denominator > 0;
Rs(valid) = currentIncidence(valid) ./ denominator(valid);

% Share of the generation-interval distribution inside the curve at each time.
coverage = cumsum(w(1:n-1));
coverage = min(coverage,1);

if minCoverage>0
    Rs(coverage < minCoverage - 1e-9) = NaN;
end
end
