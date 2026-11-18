function [PII, rangeFrac, scale, PIIbox, ratio] = subepidemicPII(lo, hi, LBbox, UBbox, scaleUser)
%SUBEPIDEMICPII  Parameter identifiability index (PII) and its companions.
%
%   [PII, rangeFrac, scale, PIIbox, ratio] = subepidemicPII(lo, hi, LBbox, UBbox)
%   [...] = subepidemicPII(lo, hi, LBbox, UBbox, scaleUser)
%
%   The same definitions as quantdiffPII.m of the QuantDiffForecast toolbox,
%   so that the two toolboxes report the same quantities.
%
%   lo, hi        2.5% and 97.5% bootstrap quantiles of each parameter
%   LBbox, UBbox  search bounds of the same parameters (subepidemicSearchSpace.m)
%   scaleUser     optional override of SCALE per parameter: 1 = log scale,
%                 0 = range fraction, NaN = automatic. Shorter vectors are
%                 padded with NaN.
%
%   PII        log10(hi/lo), the width of the 95% interval in decades (1 = the
%              interval spans a factor of 10). Defined only when lo > 0;
%              NaN otherwise.
%   rangeFrac  (hi - lo)/(UBbox - LBbox): the share of the search range that
%              the interval still covers, in [0,1]. Defined for every free
%              parameter; NaN for a fixed one (UBbox == LBbox). 1 = the data
%              add nothing to the bounds. It depends on how wide the bounds
%              were set, so read it together with PII, not instead of it.
%   scale      which measure is the headline for each parameter:
%                1  log scale (PII): the lower bound is positive, so relative
%                   precision is meaningful (r, K, alpha, d)
%                0  range fraction: the lower bound is <= 0 (p, and a of the
%                   GRM and Richards kernels), where a log ratio is undefined
%                   or inflated
%   PIIbox     log10(UBbox/LBbox), the largest PII the bounds allow; NaN when
%              LBbox <= 0.
%   ratio      bound-censoring diagnostic, ~1 = the interval is the search range:
%                scale 1:  PII ./ PIIbox   (NaN for fixed parameters)
%                scale 0:  rangeFrac
%
%   See also subepidemicIdentifiability, subepidemicBoundHitFractions,
%   subepidemicQuantileMCSE, sci_interval
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

lo = double(full(lo(:)')); hi = double(full(hi(:)'));
LBbox = double(full(LBbox(:)')); UBbox = double(full(UBbox(:)'));
n = numel(lo);
if numel(hi) ~= n || numel(LBbox) ~= n || numel(UBbox) ~= n
    error('SubEpiPredict:PIISize', 'lo, hi, LBbox and UBbox must have the same number of elements.');
end
if ~isreal(lo) || ~isreal(hi) || ~isreal(LBbox) || ~isreal(UBbox)
    error('SubEpiPredict:PIIComplex', 'Inputs must be real.');
end

% Log-scale index: only where both quantiles are positive.
PII = nan(1, n);
ok = lo > 0 & hi > 0;
PII(ok) = log10(hi(ok)) - log10(lo(ok));

% Share of the search range covered by the interval.
width = UBbox - LBbox;
rangeFrac = nan(1, n);
free = width > 0 & isfinite(width) & isfinite(lo) & isfinite(hi);
rangeFrac(free) = min(max((hi(free) - lo(free)) ./ width(free), 0), 1);

% Headline scale: automatic rule, then the user's override where given.
scale = double(LBbox > 0);
if nargin >= 5 && ~isempty(scaleUser)
    scaleUser = double(scaleUser(:)');
    m = min(numel(scaleUser), n);
    given = false(1, n);
    given(1:m) = isfinite(scaleUser(1:m));
    if any(scaleUser(given(1:m)) ~= 0 & scaleUser(given(1:m)) ~= 1)
        error('SubEpiPredict:PIIScale', 'PIIscale entries must be 1 (log scale), 0 (range fraction) or NaN (automatic).');
    end
    scale(given) = scaleUser(given);
end

% Largest PII the bounds allow; undefined for a nonpositive lower bound.
PIIbox = nan(1, n);
pos = LBbox > 0 & UBbox > 0;
PIIbox(pos) = log10(UBbox(pos)) - log10(LBbox(pos));

% Bound-censoring diagnostic on the headline scale of each parameter.
ratio = nan(1, n);
isLog = scale == 1;
ratio(isLog)  = PII(isLog) ./ PIIbox(isLog);
ratio(~isLog) = rangeFrac(~isLog);
ratio(isLog & PIIbox == 0) = NaN;          % fixed parameter: 0/0, made explicit
end
