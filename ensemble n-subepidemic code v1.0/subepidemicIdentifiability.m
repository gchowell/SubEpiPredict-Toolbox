function pii = subepidemicIdentifiability(Phatss, LB, UB, npatches, P0, scaleByType, B2, boundHitTol)
%SUBEPIDEMICIDENTIFIABILITY  PII metrics of one ranked sub-epidemic model.
%
%   pii = subepidemicIdentifiability(Phatss, LB, UB, npatches, P0)
%   pii = subepidemicIdentifiability(Phatss, LB, UB, npatches, P0, scaleByType, B2, boundHitTol)
%
%   The identifiability measures that the QuantDiffForecast toolbox reports
%   for each parameter (quantdiffPII.m, boundHitFractions.m, quantileMCSE),
%   computed from the bootstrap refits of one ranked model.
%
%   Inputs
%     Phatss       M-by-(4n+2) bootstrap refits, columns
%                  [r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]
%     LB, UB       search bounds in the same layout (subepidemicSearchSpace.m)
%     npatches     n, the number of sub-epidemics of the model
%     P0           the ranked fit to the observed data, same layout ([] = NaN)
%     scaleByType  optional headline scale per parameter type, [r p a K alpha d]:
%                  1 = log scale, 0 = range fraction, NaN = automatic
%                  (setting PIIscale of options_fit_settings.m). [] = automatic.
%     B2           resamples of the Monte Carlo error (default 500)
%     boundHitTol  relative tolerance of the bound-hit fractions (default 1e-3)
%
%   Output: a struct with one entry per reported parameter (k-by-1 fields).
%   r, p, a and K of every sub-epidemic are always listed; alpha and d only
%   when the estimation method estimates them (LB < UB).
%     name           'r_1', 'p_2', 'K_1', 'alpha', 'd', ...
%     type           'r', 'p', 'a', 'K', 'alpha' or 'd'
%     subepidemic    index of the sub-epidemic (NaN for alpha and d)
%     column         column of Phatss
%     estimate       value in the ranked fit
%     median         median of the refits
%     CI_LB, CI_UB   2.5% and 97.5% quantiles of the refits
%     PII            log10(CI_UB/CI_LB), interval width in decades; NaN when CI_LB <= 0
%     PII_rangeFrac  (CI_UB - CI_LB)/(UB - LB), in [0,1]; NaN for a fixed parameter
%     PII_scale      headline measure: 1 = PII (LB > 0), 0 = PII_rangeFrac (LB <= 0)
%     PII_box        log10(UB/LB), the largest PII the bounds allow; NaN when LB <= 0
%     PII_ratio      ~1 = the interval is the search range: PII/PII_box on the
%                    log scale, PII_rangeFrac otherwise; NaN for a fixed parameter
%     fracLB, fracUB share of the refits within boundHitTol (relative) of LB, UB
%     MCSE_CI_LB, MCSE_CI_UB, MCSE_PII
%                    Monte Carlo SE of the interval limits and of the PII
%     bound_LB, bound_UB  the search bounds
%     fixed          true where LB == UB (the parameter is not estimated)
%   and the scalars M (refits), B2 and boundHitTol.
%
%   Reading the numbers: PII near 0 is a tight interval, PII above 1 an
%   interval wider than a factor of 10. An interval is evidence about the
%   data only when PII_ratio is well below 1 and fracLB and fracUB are near
%   0; otherwise the search bounds, not the data, set its limits. A
%   difference between two PII values smaller than about twice MCSE_PII is
%   bootstrap noise.
%
%   See also subepidemicPII, subepidemicBoundHitFractions,
%   subepidemicQuantileMCSE, subepidemicPrintPII, subepidemicWritePII
%
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SubEpiPredict:InvalidPatchCount','npatches must be a positive integer.');
end
n = full(double(npatches));
nz = 4*n+2;
Phatss = double(full(Phatss));
LB = double(full(LB(:)'));
UB = double(full(UB(:)'));
if size(Phatss,2)~=nz || numel(LB)~=nz || numel(UB)~=nz
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Phatss, LB and UB must have %d columns for %d sub-epidemics.',nz,n);
end
if nargin<5 || isempty(P0)
    P0 = nan(1,nz);
end
P0 = double(full(P0(:)'));
if numel(P0)~=nz
    error('SubEpiPredict:ParameterDimensionMismatch','P0 must have %d entries.',nz);
end
if nargin<6 || isempty(scaleByType)
    scaleByType = nan(1,6);
end
scaleByType = double(scaleByType(:)');
if numel(scaleByType)~=6
    error('SubEpiPredict:PIIScale', ...
        'PIIscale must have 6 entries, one per parameter type: [r p a K alpha d].');
end
if nargin<7 || isempty(B2)
    B2 = 500;
end
if nargin<8 || isempty(boundHitTol)
    boundHitTol = 1e-3;
end

% One scale entry per column of Phatss.
scaleUser = [repmat(scaleByType(1),1,n) repmat(scaleByType(2),1,n) ...
    repmat(scaleByType(3),1,n) repmat(scaleByType(4),1,n) scaleByType(5) scaleByType(6)];

lo = quantile(Phatss,0.025,1);
hi = quantile(Phatss,0.975,1);
med = median(Phatss,1);

[PII,rangeFrac,scale,PIIbox,ratio] = subepidemicPII(lo,hi,LB,UB,scaleUser);
[fracLB,fracUB] = subepidemicBoundHitFractions(Phatss,LB,UB,boundHitTol);
[mcseLB,mcseUB,mcsePII] = subepidemicQuantileMCSE(Phatss,B2);

% Parameters reported: r, p, a, K of every sub-epidemic; alpha and d when estimated.
keep = [1:4*n, 4*n+find(UB(4*n+1:end)>LB(4*n+1:end))];

letters = 'rpaK';
names = cell(nz,1); types = cell(nz,1); sub = nan(nz,1);
for j = 1:4*n
    block = ceil(j/n);
    sub(j) = j-(block-1)*n;
    types{j} = letters(block);
    names{j} = sprintf('%s_%d',letters(block),sub(j));
end
names{4*n+1} = 'alpha'; types{4*n+1} = 'alpha';
names{4*n+2} = 'd';     types{4*n+2} = 'd';

pii = struct();
pii.name = names(keep);
pii.type = types(keep);
pii.subepidemic = sub(keep);
pii.column = keep(:);
pii.estimate = P0(keep)';
pii.median = med(keep)';
pii.CI_LB = lo(keep)';
pii.CI_UB = hi(keep)';
pii.PII = PII(keep)';
pii.PII_rangeFrac = rangeFrac(keep)';
pii.PII_scale = scale(keep)';
pii.PII_box = PIIbox(keep)';
pii.PII_ratio = ratio(keep)';
pii.fracLB = fracLB(keep)';
pii.fracUB = fracUB(keep)';
pii.MCSE_CI_LB = mcseLB(keep)';
pii.MCSE_CI_UB = mcseUB(keep)';
pii.MCSE_PII = mcsePII(keep)';
pii.bound_LB = LB(keep)';
pii.bound_UB = UB(keep)';
pii.fixed = (LB(keep)==UB(keep))';
pii.M = size(Phatss,1);
pii.B2 = B2;
pii.boundHitTol = boundHitTol;
end
