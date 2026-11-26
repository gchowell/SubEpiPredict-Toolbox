function [mcseLB, mcseUB, mcsePII] = subepidemicQuantileMCSE(Phatss, B2)
%SUBEPIDEMICQUANTILEMCSE  Monte Carlo standard error of the bootstrap interval and PII.
%
%   [mcseLB, mcseUB, mcsePII] = subepidemicQuantileMCSE(Phatss, B2)
%
%   The bootstrap-the-bootstrap of the QuantDiffForecast toolbox (quantileMCSE
%   in Run_Fit_ODEModel.m). The rows of Phatss (M-by-d bootstrap estimates)
%   are resampled with replacement B2 times (default 500); the 2.5% and 97.5%
%   quantiles and PII = log10(q97.5/q2.5) of every column are recomputed, and
%   the SD across resamples is returned: the seed-to-seed noise of the
%   interval limits and of the PII that comes from using M refits and not
%   infinitely many. No refits are performed.
%
%   Columns with zero spread (fixed parameters) and all-NaN columns give NaN.
%   mcsePII is NaN where fewer than two resamples have positive quantiles
%   (log10 undefined), as for the PII itself.
%
%   The resampling uses its own random stream with a fixed seed, so the
%   result is repeatable and the global generator (setting rngSeed) is not
%   advanced.
%
%   See also subepidemicPII, subepidemicIdentifiability
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if nargin < 2 || isempty(B2)
    B2 = 500;
end
Phatss = double(full(Phatss));
[M, d] = size(Phatss);
mcseLB = nan(1, d); mcseUB = nan(1, d); mcsePII = nan(1, d);
if M < 2 || B2 < 2, return, end
rs  = RandStream('mt19937ar', 'Seed', 20261030);           % local stream: the global RNG is untouched
idx = randi(rs, M, M, B2);                                 % one resample per column
for k = 1:d
    x = Phatss(:, k);
    if all(isnan(x)) || max(x) - min(x) == 0, continue, end   % fixed / degenerate
    xr  = reshape(x(idx), M, B2);                          % M-by-B2 resampled estimates
    q   = quantile(xr, [0.025 0.975], 1);                  % 2-by-B2
    mcseLB(k) = std(q(1, :), 0, 2);
    mcseUB(k) = std(q(2, :), 0, 2);
    ok = q(1, :) > 0 & q(2, :) > 0;
    if nnz(ok) >= 2
        mcsePII(k) = std(log10(q(2, ok)) - log10(q(1, ok)), 0, 2);
    end
end
end
