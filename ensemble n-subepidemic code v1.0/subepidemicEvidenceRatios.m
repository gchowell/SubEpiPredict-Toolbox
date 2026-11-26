function [evidenceRatio, weak] = subepidemicEvidenceRatios(AICc, topmodels, area, limit)
%SUBEPIDEMICEVIDENCERATIOS Evidence ratios of the top-ranked models, with a warning for weak ones.
%   [evidenceRatio, weak] = subepidemicEvidenceRatios(AICc, topmodels, area, limit)
%
%   AICc       AICc of the ranked candidates, best first (column 3 of the
%              ranking of the screen)
%   topmodels  number of top-ranked models kept (topmodelsx)
%   area       the area (column of the data file), for the message
%   limit      evidence ratio from which a model counts as weakly supported
%              (setting evidenceRatioLimit, default 100; Inf: no warning)
%
%   evidenceRatio(k) = exp((AICc(k) - AICc(1))/2), k = 1..topmodels: how many
%   times more likely ranked model 1 is than ranked model k, given the data.
%   A ratio of 100 corresponds to an AICc difference of 9.2. weak lists the
%   ranked models with a ratio of at least LIMIT.
%
%   REPORT ONLY. Weak models are not dropped: the number of top models stays
%   topmodelsx, so the ranked files and the ensembles are the same as without
%   this check. The warning SubEpiPredict:WeakTopModel says which models the
%   data hardly support. It matters most for ensembles with equal weights
%   (weight_type1 = -1), where a weak model counts as much as the best one;
%   Akaike weights (weight_type1 = 1) already give it almost nothing.
%
%   Added in revision l (2026-11-18). It replaces a line of
%   Run_Fit_subepidemicFramework.m that computed the models with a ratio
%   below 100 and did not use the result.

if nargin<4 || isempty(limit)
    limit = 100;
end
% full(): the ranking of the screen is a sparse matrix, and MATLAB's arrayfun
% (used for the message below) refuses sparse input.
AICc = full(double(AICc(:)));
topmodels = min(max(0,floor(full(double(topmodels)))), numel(AICc));
evidenceRatio = nan(topmodels,1);
weak = zeros(0,1);
if topmodels==0
    return
end
if ~isfinite(AICc(1))
    % no usable best model: there is nothing to compare with
    return
end
delta = AICc(1:topmodels) - AICc(1);
evidenceRatio = exp(delta/2);               % Inf where AICc is Inf (too few data points for the model)
evidenceRatio(1) = 1;
weak = find(evidenceRatio >= limit);
if isempty(weak)
    return
end

ratioText = strjoin(arrayfun(@localNumber, evidenceRatio(weak)', 'UniformOutput', false), ', ');
deltaText = strjoin(arrayfun(@(v) sprintf('%.1f', v), delta(weak)', 'UniformOutput', false), ', ');
if numel(weak)==1
    who = sprintf('ranked model %d has an evidence ratio', weak);
    them = 'it';
else
    who = sprintf('ranked models %s have evidence ratios', strjoin(arrayfun(@num2str, weak', 'UniformOutput', false), ', '));
    them = 'them';
end
subepidemicWarning('SubEpiPredict:WeakTopModel', ...
    ['Area %d: %s of %s against ranked model 1 (AICc higher by %s). The data give %s little support ', ...
     '(evidence ratio of %s or more). %s kept, because topmodelsx = %d; an ensemble with equal weights ', ...
     '(weight_type1 = -1) counts %s as much as the best model.'], ...
    area, who, ratioText, deltaText, them, localNumber(limit), localCapital(them), topmodels, them);
end

% -----------------------------------------------------------------------------
function s = localNumber(v)
if isinf(v)
    s = 'Inf';
elseif v < 1e4
    s = sprintf('%.3g', v);
else
    s = sprintf('%.2g', v);
end
end

function s = localCapital(them)
if strcmp(them,'it')
    s = 'It is';
else
    s = 'They are';
end
end
