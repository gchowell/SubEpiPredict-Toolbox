function [LB,UB,isLog,z0] = subepidemicSearchSpace(npatches,data,flag1,method1,settings)
%SUBEPIDEMICSEARCHSPACE Bounds, coordinate scales and seed of one candidate.
%   [LB,UB,isLog,z0] = subepidemicSearchSpace(npatches,data,flag1,method1,settings)
%
%   Layout of every vector: [r(1:n), p(1:n), a(1:n), K(1:n), alpha, d],
%   4*n+2 entries, as in unpackSubepidemicParameters.
%
%   LB, UB  With default settings these are the bounds the toolbox used
%           before (getbounds.m for r, p, a and K, and the alpha and d
%           bounds of the fitting drivers), with one exception: the upper
%           bound of r is never below settings.rUpperBoundFactor (default
%           10) times the largest count. rUpperBoundFactor = 0 gives the
%           legacy bounds exactly, apart from the Gompertz kernel (below).
%           A parameter with LB == UB is fixed. The
%           fields rLB ... KUB of SETTINGS replace a bound (see
%           subepidemicFitSettings.m). Reversed bounds are an error.
%           Gompertz kernel (flag1 = 5): a is fixed at 1. getbounds.m let
%           it vary from 0 to 10, but the Gompertz equation does not contain
%           a, so the search had a dimension with no effect on the fit, and
%           the histogram and interval of a were noise. AICc already counted
%           two parameters per Gompertz sub-epidemic (get_nparams.m).
%   isLog   true for the coordinates searched on a log scale: r, K and
%           alpha, when the setting is on, the lower bound is positive and
%           the parameter is free.
%   z0      the legacy seed (r = 0.1, p = 0.9 or 1, a = 1, K = sum of the
%           data, alpha = 1, d = 1), projected into the bounds.
%
%   DATA is the calibration incidence as a column (the unsmoothed series the
%   drivers pass to getbounds). FLAG1 is the growth-kernel code 0 to 5 and
%   METHOD1 the estimation method 0 to 6. No global variables are used.
%   Method 6 (sum of absolute deviations) gets fixed alpha and d, like
%   methods 0 to 2; the legacy switch had no case for it.
%
%   See also subepidemicFitSettings, fitSubepidemicModel, getbounds

if nargin<5 || isempty(settings)
    settings = subepidemicFitSettings();
end
if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SubEpiPredict:InvalidPatchCount','npatches must be a positive integer.');
end
npatches = full(double(npatches));
if ~(isnumeric(data) && isreal(data) && isvector(data) && numel(data)>=2 && all(isfinite(data(:))))
    error('SubEpiPredict:InvalidData', ...
        'data must be a finite real vector with at least two observations.');
end
data = double(data(:));
if ~(isnumeric(flag1) && isscalar(flag1) && any(flag1==0:5))
    error('SubEpiPredict:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
end
if ~(isnumeric(method1) && isscalar(method1) && any(method1==0:6))
    error('SubEpiPredict:InvalidMethod','method1 must be an integer from 0 to 6.');
end

% ---- legacy bounds of getbounds.m (same numbers, no globals) ----------------
rlb = mean(abs(data(1:2)))/200;
rub = max(abs(data(1:2)))*50;
% The legacy upper bound of r depends on the first two observations only: 50
% for a series that starts with 1, 0. Late-onset candidates and bootstrap
% refits reached it, which cut the bootstrap interval of r off at the bound.
% A sub-epidemic that starts at a cumulative count of 1 or more has an
% incidence of about r or more at its start, so r far above the largest
% count cannot fit count data; a multiple of the largest count is therefore a
% bound that scales with the data and does not constrain a sensible fit.
if isfield(settings,'rUpperBoundFactor') && ~isempty(settings.rUpperBoundFactor)
    rub = max(rub,settings.rUpperBoundFactor*max(abs(data)));
end
Kmax = 10000000000;
switch flag1
    case 0 % GGM
        pb = [0 1]; ab = [1 1]; Kb = [20 20];
    case 1 % GLM
        pb = [0 1]; ab = [1 1]; Kb = [20 Kmax];
    case 2 % GRM
        pb = [0 1]; ab = [0 10]; Kb = [20 Kmax];
    case 3 % Logistic
        pb = [1 1]; ab = [1 1]; Kb = [20 Kmax];
    case 4 % Richards
        pb = [1 1]; ab = [0 10]; Kb = [20 Kmax];
    case 5 % Gompertz: dC/dt = r*C*log(K/C) has no a (see modifiedLogisticGrowthPatch.m)
        pb = [1 1]; ab = [1 1]; Kb = [20 Kmax];
end
rLB = localBound(settings,'rLB',rlb,npatches);
rUB = localBound(settings,'rUB',rub,npatches);
pLB = localBound(settings,'pLB',pb(1),npatches);
pUB = localBound(settings,'pUB',pb(2),npatches);
aLB = localBound(settings,'aLB',ab(1),npatches);
aUB = localBound(settings,'aUB',ab(2),npatches);
KLB = localBound(settings,'KLB',Kb(1),npatches);
KUB = localBound(settings,'KUB',Kb(2),npatches);

% ---- legacy bounds of the error-structure parameters (alpha, d) -------------
switch method1
    case {0,1,2,6}
        LBe = [0 0];      UBe = [0 0];
    case {3,4}
        LBe = [10^-8 0];  UBe = [10^5 0];
    case 5
        LBe = [10^-8 0.2]; UBe = [10^5 10^2];
end

LB = [rLB pLB aLB KLB LBe];
UB = [rUB pUB aUB KUB UBe];
if any(LB>UB)
    error('SubEpiPredict:ReversedBounds', ...
        'A lower bound exceeds its upper bound. Check the bound settings (rLB ... KUB).');
end

% ---- coordinate scales ------------------------------------------------------
isLog = false(1,4*npatches+2);
isLog(1:npatches) = settings.logScaleR;
isLog(3*npatches+1:4*npatches) = settings.logScaleK;
isLog(4*npatches+1) = settings.logScaleAlpha;
isLog = isLog & LB>0 & UB>LB;

% ---- legacy seed ------------------------------------------------------------
I0 = data(1);
total = sum(data);
rs = zeros(1,npatches)+0.1;
ps = zeros(1,npatches)+0.9;
as = ones(1,npatches);
Ks = zeros(1,npatches)+total;
if flag1==3 || flag1==4 || flag1==5 % logistic, Richards and Gompertz have p = 1
    ps(:) = 1;
end
if flag1==5
    rs(:) = 1-I0/total;
    as(:) = rs(1)/log(total/I0);
end
z0 = [rs ps as Ks 1 1];
z0 = min(max(z0,LB),UB);
z0(~isfinite(z0)) = LB(~isfinite(z0));
end

% -----------------------------------------------------------------------------
function b = localBound(settings,name,legacy,npatches)
% One bound per sub-epidemic: the legacy value, or the setting when given.
b = zeros(1,npatches)+legacy;
if isfield(settings,name) && ~isempty(settings.(name))
    v = double(settings.(name)(:).');
    if isscalar(v)
        b(:) = v;
    elseif numel(v)>=npatches
        b = v(1:npatches);
    else
        error('SubEpiPredict:FitSettings', ...
            'Setting %s has %d entries; %d sub-epidemics need a scalar or at least %d values.', ...
            name,numel(v),npatches,npatches);
    end
end
end
