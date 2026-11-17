function [value,detail] = subepidemicObjective(z,context)
%SUBEPIDEMICOBJECTIVE Fitting objective with an explicit context, no globals.
%   [value,detail] = subepidemicObjective(z,context)
%
%   z        parameter vector [r(1:n), p(1:n), a(1:n), K(1:n), alpha, d]
%   context  struct with fields
%              timevect, ydata   observation times and the series fitted
%              I0                initial value (first observation)
%              npatches          number of sub-epidemics n
%              onset_fixed       0 asynchronous, 1 synchronous
%              onset_thr         onset threshold C_thr
%              flag1             growth-kernel code 0 to 5
%              method1           0 least squares, 1 Poisson, 3/4/5 negative
%                                binomial, 6 sum of absolute deviations
%              odeOptions        odeset structure (subepidemicOdeOptions)
%
%   value    the same number the legacy objective
%            plotModifiedLogisticGrowthPatchMethodLogLik returns for a valid
%            trajectory: SSE, the negative log-likelihood without its
%            data-only constant, or the sum of absolute deviations. AICc is
%            computed from it as before.
%   detail   struct: success, message, nActive (sub-epidemics active by the
%            end of the grid), activation (fifth output of
%            simulateSubepidemic), curve (fitted incidence).
%
%   A parameter vector is INVALID (detail.success = false, value = 1e10, the
%   legacy penalty) when the simulation fails, the trajectory is identically
%   zero, or the objective is not a finite real number. Unlike the legacy
%   wrapper, a valid trajectory with a very large SSE is not penalized, so
%   the optimizer keeps a slope to follow from a poor start point.
%
%   The legacy wrappers are unchanged and still read their globals; the
%   fitting drivers no longer call them.
%
%   See also fitSubepidemicModel, simulateSubepidemic, subepidemicNegativeBinomialNLL

value = 1e10;
detail = struct('success',false,'message','','nActive',0,'activation',[],'curve',[]);
try
    n = context.npatches;
    if ~(isnumeric(z) && isreal(z) && numel(z)==4*n+2 && all(isfinite(z(:))))
        detail.message = 'Parameter vector has the wrong size or nonfinite values.';
        return
    end
    [rs1,ps1,as1,Ks1,alpha,d] = unpackSubepidemicParameters(z,n);
    [~,~,yfit,~,activation] = simulateSubepidemic(context.timevect,context.I0,n, ...
        context.onset_fixed,context.onset_thr,context.flag1,rs1,ps1,as1,Ks1, ...
        context.odeOptions);
    detail.activation = activation;
    detail.nActive = sum(activation.active);
    if sum(yfit)==0
        detail.message = 'The fitted trajectory is identically zero.';
        return
    end
    yfit(yfit==0) = 0.001; % legacy convention, kept so objectives are comparable
    ydata = context.ydata(:);
    yfit = yfit(:);
    switch context.method1
        case 0
            v = sum((ydata-yfit).^2);
        case 1
            v = -sum(ydata.*log(yfit)-yfit);
        case {3,4,5}
            v = subepidemicNegativeBinomialNLL(ydata,yfit,alpha,d,context.method1);
        case 6
            v = sum(abs(ydata-yfit));
        otherwise
            error('SubEpiPredict:InvalidMethod','method1 must be 0, 1, 3, 4, 5 or 6.');
    end
    if ~(isscalar(v) && isreal(v) && isfinite(v))
        detail.message = 'The objective is not a finite real number at this parameter vector.';
        return
    end
    value = v;
    detail.success = true;
    detail.curve = yfit;
catch ME
    value = 1e10;
    detail.success = false;
    detail.message = sprintf('%s: %s',ME.identifier,ME.message);
end
end
