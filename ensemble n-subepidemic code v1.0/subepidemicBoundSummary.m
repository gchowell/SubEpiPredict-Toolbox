function [report,lines,nRowsAtBound] = subepidemicBoundSummary(P,LB,UB,isLog,npatches)
%SUBEPIDEMICBOUNDSUMMARY Parameters that sit at a bound, and in how many fits.
%   [report,lines] = subepidemicBoundSummary(P,LB,UB,isLog,npatches)
%
%   P holds one fitted parameter vector per ROW, in the layout
%   [r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]; LB, UB and isLog are the
%   outputs of subepidemicSearchSpace for the same model. A parameter is "at
%   a bound" by the rule of fitSubepidemicModel: its unit coordinate is
%   within 1e-6 of 0 or 1 (for r, K and alpha that is a distance on the log
%   scale). Fixed parameters (LB == UB) are never reported. Rows with a
%   non-finite entry are ignored.
%
%   report  struct array, one element per parameter and side with at least
%           one fit at the bound, in parameter order:
%             name           'r1', 'p2', 'K1', 'alpha', 'd', ...
%             index          column of P
%             side           'lower' or 'upper'
%             limit          the bound
%             count, total   fits at the bound, fits examined
%             fraction       count/total
%             setting        the field of options_fit_settings.m that moves
%                            the bound ('rUB', 'KLB', ...; '' for alpha, d)
%             isSearchLimit  true for r and K. Their bounds only limit the
%                            search, so an estimate or an interval that ends
%                            there says nothing about the data. The bounds of
%                            p and a belong to the growth model itself.
%   lines   cell array of char, one printable line per element, such as
%           'r2 at its upper limit (50) in 6 of 40'.
%
%   A percentile interval of a parameter ends exactly at a bound when more
%   than 2.5% of the bootstrap refits are there.
%
%   See also subepidemicSearchSpace, fitSubepidemicModel, subepidemicUnitMap

if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SubEpiPredict:InvalidPatchCount','npatches must be a positive integer.');
end
npatches = full(double(npatches));
nz = 4*npatches+2;
LB = double(full(LB(:).'));
UB = double(full(UB(:).'));
isLog = logical(full(isLog(:).'));
if numel(LB)~=nz || numel(UB)~=nz || numel(isLog)~=nz
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Bounds and scale flags must have %d entries for %d sub-epidemics.',nz,npatches);
end
P = double(full(P));
if isvector(P)
    P = P(:).';
end
if size(P,2)~=nz
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Each row must hold %d parameters; received %d columns.',nz,size(P,2));
end
P = P(all(isfinite(P),2),:);
total = size(P,1);

report = struct('name',{},'index',{},'side',{},'limit',{},'count',{},'total',{}, ...
    'fraction',{},'setting',{},'isSearchLimit',{});
lines = {};
nRowsAtBound = 0;
if total==0
    return
end

map = subepidemicUnitMap('make',struct('z0',LB,'LB',LB,'UB',UB,'isLog',isLog));
% Revision m: the rule of the PII table (subepidemicBoundHitFractions.m), a
% value within 0.1% of the bound, so the bootstrap summary and the columns
% atLB/atUB of the PII table give the same counts. Until revision l this
% function counted values within 1e-6 of the search range.
[~,~,hitLB,hitUB] = subepidemicBoundHitFractions(P,LB,UB);
nRowsAtBound = sum(any(hitLB(:,map.free) | hitUB(:,map.free),2));

letters = 'rpaK';
for f = 1:map.nFree
    j = map.free(f);
    if j<=4*npatches
        block = ceil(j/npatches);
        name = sprintf('%s%d',letters(block),j-(block-1)*npatches);
        prefix = letters(block);
        searchLimit = block==1 || block==4;
    elseif j==4*npatches+1
        name = 'alpha'; prefix = ''; searchLimit = false;
    else
        name = 'd'; prefix = ''; searchLimit = false;
    end
    counts = [sum(hitLB(:,j)) sum(hitUB(:,j))];
    sides = {'lower','upper'};
    limits = [LB(j) UB(j)];
    suffix = {'LB','UB'};
    for s = 1:2
        if counts(s)>0
            if isempty(prefix)
                setting = '';
            else
                setting = [prefix suffix{s}];
            end
            report(end+1) = struct('name',name,'index',j,'side',sides{s},'limit',limits(s), ...
                'count',counts(s),'total',total,'fraction',counts(s)/total, ...
                'setting',setting,'isSearchLimit',searchLimit); %#ok<AGROW>
            lines{end+1} = sprintf('%s at its %s limit (%s) in %d of %d', ...
                name,sides{s},num2str(limits(s),6),counts(s),total); %#ok<AGROW>
        end
    end
end
end
