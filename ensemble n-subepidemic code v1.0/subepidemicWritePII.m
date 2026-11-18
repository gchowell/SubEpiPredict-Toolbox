function subepidemicWritePII(filename, pii, ranks)
%SUBEPIDEMICWRITEPII  Save the PII metrics as a csv file.
%   subepidemicWritePII(filename, pii) writes the struct returned by
%   subepidemicIdentifiability.m, one row per parameter.
%
%   subepidemicWritePII(filename, piiList, ranks) writes several ranked
%   models into one file. piiList is a cell array of such structs and ranks
%   the rank of each; the rows start with two more columns, ranked_model and
%   n_subepidemics.
%
%   Columns
%     parameter, subepidemic   'r_1', 'K_2', 'alpha', ...; the sub-epidemic
%                              index (NaN for alpha and d)
%     estimate, median         ranked fit; median of the bootstrap refits
%     CI_LB, CI_UB             95% bootstrap interval
%     PII                      log10(CI_UB/CI_LB); NaN when CI_LB <= 0
%     PII_rangeFrac            (CI_UB - CI_LB)/(bound_UB - bound_LB)
%     PII_scale                1 = PII is the headline, 0 = PII_rangeFrac
%     PII_box                  log10(bound_UB/bound_LB); NaN when bound_LB <= 0
%     PII_ratio                PII/PII_box (scale 1) or PII_rangeFrac (scale 0)
%     fracLB, fracUB           share of the refits at the lower, upper bound
%     MCSE_CI_LB, MCSE_CI_UB, MCSE_PII   Monte Carlo SE of the limits and the PII
%     bound_LB, bound_UB       search bounds
%     fixed                    1 = not estimated (bound_LB = bound_UB)
%
%   See also subepidemicIdentifiability, subepidemicPrintPII

if isstruct(pii)
    pii = {pii};
end
withRank = nargin>=3 && ~isempty(ranks);
if withRank && numel(ranks)~=numel(pii)
    error('SubEpiPredict:PIIRanks','ranks must have one entry per PII struct.');
end

cols = {'estimate','median','CI_LB','CI_UB','PII','PII_rangeFrac','PII_scale','PII_box', ...
    'PII_ratio','fracLB','fracUB','MCSE_CI_LB','MCSE_CI_UB','MCSE_PII','bound_LB','bound_UB'};

fid = fopen(filename,'w');
if fid<0
    error('SubEpiPredict:PIIFile','The file %s could not be opened for writing.',filename);
end
closeFile = onCleanup(@() fclose(fid));

if withRank
    fprintf(fid,'ranked_model,n_subepidemics,');
end
fprintf(fid,'parameter,subepidemic');
fprintf(fid,',%s',cols{:});
fprintf(fid,',fixed\n');

for m = 1:numel(pii)
    S = pii{m};
    nSub = max(S.subepidemic(~isnan(S.subepidemic)));
    for k = 1:numel(S.name)
        if withRank
            fprintf(fid,'%d,%d,',ranks(m),nSub);
        end
        fprintf(fid,'%s,%s',S.name{k},localNumber(S.subepidemic(k)));
        for c = 1:numel(cols)
            fprintf(fid,',%s',localNumber(S.(cols{c})(k)));
        end
        fprintf(fid,',%d\n',S.fixed(k));
    end
end
end

% -----------------------------------------------------------------------------
function s = localNumber(x)
% Full precision; NaN and Inf as MATLAB reads them back.
if isnan(x)
    s = 'NaN';
elseif isinf(x)
    if x>0, s = 'Inf'; else, s = '-Inf'; end
else
    s = sprintf('%.10g',x);
end
end
