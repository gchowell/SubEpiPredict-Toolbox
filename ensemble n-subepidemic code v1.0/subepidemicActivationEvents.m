function [value,isterminal,direction]=subepidemicActivationEvents(t,x,onset_thr,watchedPatches)
%SUBEPIDEMICACTIVATIONEVENTS Pure upward threshold-crossing event function.
%   WATCHEDPATCHES contains the dormant patch indices, all >=2, at the start
%   of a segment. Event k means that the predecessor of WATCHEDPATCHES(k)
%   reaches ONSET_THR. No activation flags are changed in this callback.
%   Already-satisfied thresholds are handled outside the solver, before a
%   segment starts. Event indices map to WATCHEDPATCHES, not directly to x.

value=x(watchedPatches(:)-1)-onset_thr;
isterminal=ones(numel(watchedPatches),1);
direction=ones(numel(watchedPatches),1); % upward crossings only
end
