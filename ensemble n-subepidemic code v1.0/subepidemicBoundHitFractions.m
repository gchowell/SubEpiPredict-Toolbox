function [fracLB, fracUB, hitLB, hitUB] = subepidemicBoundHitFractions(Phatss, LB, UB, relTol)
%SUBEPIDEMICBOUNDHITFRACTIONS  Share of bootstrap refits sitting on the search bounds.
%
%   [fracLB, fracUB, hitLB, hitUB] = subepidemicBoundHitFractions(Phatss, LB, UB, relTol)
%
%   The bound-censoring diagnostic of the QuantDiffForecast toolbox
%   (boundHitFractions.m), which complements PII_box / PII_ratio. The 2.5% and
%   97.5% quantiles behind the PII can hide a pile-up of estimates at ONE
%   bound (80% of the refits at LB with the rest interior gives a
%   moderate-looking PII_ratio although the estimate is set by the bound).
%   This function reports, for every column of Phatss, the fraction of
%   bootstrap refits within a relative tolerance of the lower and of the
%   upper bound.
%
%   Inputs
%     Phatss  M-by-d matrix of bootstrap estimates (rows = refits, columns = parameters).
%     LB, UB  Search bounds: 1-by-d vectors (one box for all refits) or M-by-d matrices.
%     relTol  Relative tolerance (default 1e-3, i.e. 0.1% of the bound value).
%             A refit P counts as hitting LB when  P <= LB + relTol*|LB|
%             and as hitting UB when               P >= UB - relTol*|UB|.
%             For a bound that is exactly 0 the tolerance is 0 (exact hit only).
%
%   Outputs
%     fracLB, fracUB  1-by-d fractions in [0,1]; NaN for fixed columns
%                     (LB == UB) and for columns whose refits are all NaN.
%     hitLB, hitUB    M-by-d logical masks of the individual hits.
%
%   Reading the output
%     fracLB+fracUB ~ 0   the bootstrap distribution is interior; the bounds
%                         are not shaping the interval.
%     fracLB or fracUB large (e.g. > 0.05) with a moderate PII_ratio: the
%                         quantile-based PII understates the censoring; the
%                         estimate is pinned on one side by the bound.
%     both large          refits split between the two bounds, typical of a
%                         flat or bimodal likelihood over the whole range.
%
%   Since revision m subepidemicBoundSummary.m (the bootstrap summary of the
%   fit) counts with this function, so both give the same numbers. Until
%   revision l it counted a refit only within 1e-6 of the search range.
%
%   See also subepidemicPII, subepidemicIdentifiability, subepidemicBoundSummary
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 4 || isempty(relTol)
    relTol = 1e-3;
end

Phatss = double(full(Phatss));
[M, d] = size(Phatss);

LB = double(full(LB)); UB = double(full(UB));

% Broadcast a common box (given as a row or column vector) to M-by-d. Only for a
% single parameter (d==1) can a column of M values be a per-refit box; it is
% then kept as is. (Revision k: before, a column vector of bounds whose length
% happened to equal M was taken as a per-refit box for any d and refused.)
if isvector(LB) && ~(d == 1 && size(LB,1) == M), LB = repmat(LB(:)', M, 1); end
if isvector(UB) && ~(d == 1 && size(UB,1) == M), UB = repmat(UB(:)', M, 1); end

if size(LB,1) ~= M || size(UB,1) ~= M || size(LB,2) < d || size(UB,2) < d
    error('SubEpiPredict:BoundHitSizeMismatch', ...
        'LB/UB must be 1-by-d (or longer) or M-by-d (or wider) to match Phatss (%d-by-%d).', M, d);
end

LB = LB(:, 1:d);
UB = UB(:, 1:d);

% Guard against flipped boxes.
flipMask = LB > UB;
if any(flipMask(:))
    tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
end

tolLB = relTol * abs(LB);
tolUB = relTol * abs(UB);

hitLB = Phatss <= LB + tolLB;
hitUB = Phatss >= UB - tolUB;

valid = ~isnan(Phatss);
hitLB = hitLB & valid;
hitUB = hitUB & valid;

nValid = sum(valid, 1);
fracLB = sum(hitLB, 1) ./ max(nValid, 1);
fracUB = sum(hitUB, 1) ./ max(nValid, 1);

% Fixed columns (LB==UB in every refit) and all-NaN columns carry no information.
pinned = all(LB == UB, 1);
fracLB(pinned | nValid == 0) = NaN;
fracUB(pinned | nValid == 0) = NaN;

end
