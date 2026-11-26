function weights = subepidemicEnsembleWeights(weight_type1,AICc,WISC,WISF)
%SUBEPIDEMICENSEMBLEWEIGHTS Weights of the member models of an ensemble.
%   weights = subepidemicEnsembleWeights(weight_type1,AICc,WISC,WISF)
%
%   AICc holds the AICc of the member models, one per model. WISC and WISF
%   hold their weighted interval scores (calibration; earlier forecasts) and
%   are read only for codes 2 and 3; pass [] otherwise. The result is a
%   column of nonnegative weights that sum to 1, one per member.
%
%   weight_type1
%     -1  equal weights
%      0  Akaike weights (see below)
%      1  Akaike weights: exp(-(AICc - min(AICc))/2), normalized
%      2  proportional to 1/WISC
%      3  proportional to 1/WISF
%
%   Code 0 used to be (1/AICc)/sum(1/AICc). That is not a weight. AICc is
%   defined up to an additive constant, and for the count likelihoods it is
%   a large negative number: four models with AICc near -1,091,960 all got
%   0.25, whatever their differences, and AICc values of mixed sign gave a
%   negative weight and a negative number of curves to sample. Only AICc
%   differences carry information, which is what Akaike weights use. Code 0
%   now returns them and warns once per call
%   (SubEpiPredict:WeightCodeZero); use code 1 to get the same weights
%   without the warning, or -1 for equal weights.
%
%   Scores of codes 2 and 3 must be finite and positive; anything else is an
%   error, because 1/score would be infinite or would change sign.
%
%   See also getensemblesubepidemicmodels, options_forecast

if ~(isnumeric(weight_type1) && isreal(weight_type1) && isscalar(weight_type1) && ...
        any(weight_type1==[-1 0 1 2 3]))
    error('SubEpiPredict:EnsembleWeights', ...
        'weight_type1 must be -1, 0, 1, 2 or 3 (see options_forecast.m).');
end
if ~(isnumeric(AICc) && isreal(AICc) && isvector(AICc) && ~isempty(AICc) && all(isfinite(AICc(:))))
    error('SubEpiPredict:EnsembleWeights', ...
        'AICc must hold one finite value per member model.');
end
AICc = full(double(AICc(:)));
k = numel(AICc);

switch weight_type1
    case -1
        weights = ones(k,1)/k;
    case {0,1}
        if weight_type1==0
            subepidemicWarning('SubEpiPredict:WeightCodeZero', ...
                ['weight_type1 = 0 used 1/AICc, which is not a weight (equal weights for count ', ...
                 'likelihoods, negative weights when AICc changes sign). It now gives Akaike ', ...
                 'weights, as weight_type1 = 1 does.']);
        end
        relativeLikelihood = exp(-(AICc-min(AICc))/2);
        weights = relativeLikelihood/sum(relativeLikelihood);
    case 2
        weights = localInverseScore(WISC,k,'WISC',2);
    case 3
        weights = localInverseScore(WISF,k,'WISF',3);
end

if ~(all(isfinite(weights)) && all(weights>=0) && abs(sum(weights)-1)<1e-12)
    error('SubEpiPredict:EnsembleWeights', ...
        'The ensemble weights are not finite, nonnegative numbers that sum to 1.');
end
end

% -----------------------------------------------------------------------------
function weights = localInverseScore(score,k,name,code)
if ~(isnumeric(score) && isreal(score) && numel(score)==k)
    error('SubEpiPredict:EnsembleWeights', ...
        'weight_type1 = %d needs one %s value per member model (%d).',code,name,k);
end
score = full(double(score(:)));
if ~(all(isfinite(score)) && all(score>0))
    error('SubEpiPredict:EnsembleWeights', ...
        'weight_type1 = %d needs finite, positive %s values; received %s.',code,name,mat2str(score.',6));
end
weights = (1./score)/sum(1./score);
end
