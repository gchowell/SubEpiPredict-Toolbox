function [totinc,patchIncidence]=subepidemicIncidence(x,onset_fixed)
%SUBEPIDEMICINCIDENCE Convert a full raw ODE trajectory to interval incidence.
%   X must begin at the calibration initial time and contain the raw states
%   returned by the ODE solver, not an already seed-corrected cumulative curve
%   or a forecast-only tail. Columns are patches; rows are observation times.
%   Later rows are interval increments (not instantaneous derivatives).
%   Only the first row is corrected for the artificial one-case seeds in
%   asynchronous patches 2..n. Synchronous models have no artificial seeds.
%   Sum(PATCHINCIDENCE,2) is TOTINC at every row, including the first.
%   With initializeSubepidemicState, TOTINC(1) equals I0 up to roundoff.

if ~(isnumeric(x) && ismatrix(x) && ~isempty(x))
    error('SubEpiPredict:InvalidTrajectory','x must be a nonempty numeric matrix.');
end
if ~((isnumeric(onset_fixed) || islogical(onset_fixed)) && ...
        isreal(onset_fixed) && isscalar(onset_fixed) && ...
        (onset_fixed==0 || onset_fixed==1))
    error('SubEpiPredict:InvalidOnsetMode','onset_fixed must be 0 (asynchronous) or 1 (synchronous).');
end

patchIncidence=[x(1,:);diff(x,1,1)];
if onset_fixed==0
    patchIncidence(1,2:end)=patchIncidence(1,2:end)-1;
end
totinc=sum(patchIncidence,2);
end
