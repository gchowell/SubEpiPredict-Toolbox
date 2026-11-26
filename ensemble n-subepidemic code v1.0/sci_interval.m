function s = sci_interval(lb, ub)
%SCI_INTERVAL log10 width of a 95% interval, the toolbox's identifiability score.
% The same function, with the same convention, as in the sub-epidemic SEIR
% toolbox: log10(UB / LB), the number of orders of magnitude the interval
% spans, whatever the size of the parameter.
%
% The n-sub-epidemic toolbox computed log10(UB/LB) with no guard, which is
% Inf when the lower limit is 0 (for example p of a sub-epidemic whose
% bootstrap refits reach p = 0) and NaN when both limits are 0. Here the lower
% limit is floored at floorRatio*UB, which changes nothing unless UB/LB
% exceeds 1/floorRatio = 1e6.
%
% Special cases:
%   LB == UB (a fixed parameter, zero-width interval)   -> 0
%   LB == 0 < UB (an interval that reaches zero)        -> log10(1/floorRatio) = 6,
%       the cap of the score. Read 6 as "not bounded away from zero".
%   NaN in either limit                                 -> NaN
% Vector inputs are allowed; lb and ub must have the same size. The limits are
% assumed nonnegative, as all toolbox parameters are.
%
% One line differs from the SEIR toolbox's copy: the last assignment. max()
% ignores NaN, so without it a NaN lower limit is replaced by the floor and
% the score is 6 instead of NaN.
floorRatio = 1e-6;
s = log10(ub ./ max(lb, floorRatio .* ub));
s(ub == lb) = 0;
s(isnan(lb) | isnan(ub)) = NaN;
end
