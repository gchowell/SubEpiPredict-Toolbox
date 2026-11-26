function subepidemicPrintPII(pii, heading)
%SUBEPIDEMICPRINTPII  Print the PII metrics of one ranked model.
%   subepidemicPrintPII(pii, heading) prints the struct returned by
%   subepidemicIdentifiability.m as a table, one line per parameter. Fixed
%   parameters (lower bound = upper bound) are listed as 'fixed'.
%
%   Columns: estimate and 95% interval; PII = log10(UB/LB) of the interval;
%   range = share of the search range the interval covers; scale = headline
%   measure (log: PII, range: range fraction); box = largest PII the bounds
%   allow; ratio = PII/box (log) or range fraction, ~1 when the interval is
%   the search range; atLB, atUB = share of the refits at a bound; MCSE = Monte
%   Carlo SE of the PII.
%
%   See also subepidemicIdentifiability, subepidemicWritePII

if nargin<2 || isempty(heading)
    heading = 'Parameter identifiability (PII)';
end
fprintf('\n%s (%d bootstrap refits):\n',heading,pii.M);
fprintf('  %-7s %10s %10s %10s %7s %7s %6s %7s %7s %6s %6s %7s\n', ...
    'param','estimate','CI LB','CI UB','PII','range','scale','box','ratio','atLB','atUB','MCSE');
for k = 1:numel(pii.name)
    if pii.fixed(k)
        fprintf('  %-7s %10.4g %s\n',pii.name{k},pii.estimate(k),'  fixed');
        continue
    end
    if pii.PII_scale(k)==1
        scaleName = 'log';
    else
        scaleName = 'range';
    end
    fprintf('  %-7s %10.4g %10.4g %10.4g %7.3f %7.3f %6s %7.3f %7.3f %6.3f %6.3f %7.3f\n', ...
        pii.name{k},pii.estimate(k),pii.CI_LB(k),pii.CI_UB(k),pii.PII(k),pii.PII_rangeFrac(k), ...
        scaleName,pii.PII_box(k),pii.PII_ratio(k),pii.fracLB(k),pii.fracUB(k),pii.MCSE_PII(k));
end
fprintf('  PII = log10(CI UB/CI LB); NaN when CI LB <= 0 (read the range fraction). ratio near 1 or atLB/atUB above 0.05: the bounds set the interval.\n\n');
end
