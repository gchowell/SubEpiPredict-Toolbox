function [IC,invasions0,timeinvasions0,Cinvasions0]=initializeSubepidemicState(I0,npatches,onset_fixed)
%INITIALIZESUBEPIDEMICSTATE Initial state for every subepidemic solve.
%   ONSET_FIXED, not the threshold value, defines the initialization:
%     0: asynchronous, IC=[I0;1;...;1], initially only patch 1 is active.
%     1: synchronous,  IC=repmat(I0/npatches,npatches,1), all active.
%   A zero threshold does NOT change an asynchronous model into a
%   synchronous model. simulateSubepidemic handles threshold activation
%   outside the derivative, including thresholds satisfied at the initial time.
%   The one-case seeds in later asynchronous patches are numerical initial
%   states, not additional cases in the first reported incidence observation.
%   This helper is pure: it does not read or change global variables.

if ~(isnumeric(I0) && isreal(I0) && isscalar(I0) && isfinite(I0) && I0>=0)
    error('SubEpiPredict:InvalidI0','I0 must be a finite, real, nonnegative scalar.');
end
if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SubEpiPredict:InvalidPatchCount','npatches must be a positive integer.');
end
if ~((isnumeric(onset_fixed) || islogical(onset_fixed)) && ...
        isreal(onset_fixed) && isscalar(onset_fixed) && ...
        (onset_fixed==0 || onset_fixed==1))
    error('SubEpiPredict:InvalidOnsetMode','onset_fixed must be 0 (asynchronous) or 1 (synchronous).');
end

IC=ones(npatches,1);
invasions0=zeros(npatches,1);
timeinvasions0=zeros(npatches,1);
Cinvasions0=zeros(npatches,1);

if onset_fixed==0
    IC(1)=I0;
    invasions0(1)=1;
else
    IC(:)=I0/npatches;
    invasions0(:)=1;
end
end
