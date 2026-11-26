function [ok,x,totinc,patchIncidence,message] = subepidemicSimulateRefit(timevect,I0,npatches,onset_fixed,onset_thr,flag1,rs1,ps1,as1,Ks1)
%SUBEPIDEMICSIMULATEREFIT Simulate one bootstrap refit beyond the fitted window.
%   [ok,x,totinc,patchIncidence,message] = subepidemicSimulateRefit(...) calls
%   simulateSubepidemic.m with the same arguments. ok is false, and message
%   says why, when the curve cannot be integrated to the last time or is not
%   finite. Any other error (wrong sizes, invalid settings) is not caught.
%
%   Used by the forecast and Rt scripts, which extend each refit past the
%   calibration window; an extreme parameter set can fail there although it
%   was simulated without trouble during the fit. What happens to a refit
%   that fails is decided by subepidemicSkippedRefits.m. The fit-plot and
%   ranking scripts simulate over the calibration window only, where a
%   failure would mean the saved fit does not match the code, so they still
%   stop.

x = []; totinc = []; patchIncidence = []; message = '';
try
    [~,x,totinc,patchIncidence] = simulateSubepidemic(timevect,I0,npatches,onset_fixed, ...
        onset_thr,flag1,rs1,ps1,as1,Ks1);
catch simulationError
    integrationFailures = {'SubEpiPredict:IncompleteIntegration','SubEpiPredict:InvalidTrajectory', ...
        'SubEpiPredict:InvalidSegment','SubEpiPredict:MissedActivationEvent', ...
        'SubEpiPredict:ActivationProgressFailure'};
    if ~any(strcmp(simulationError.identifier,integrationFailures))
        rethrow(simulationError)
    end
    ok = false;
    message = simulationError.message;
    return
end
ok = size(x,1)==numel(timevect) && all(isfinite(totinc(:)));
if ~ok
    message = 'the curve is incomplete or not finite';
end
end
