function [t,x,totinc,patchIncidence,activation]=simulateSubepidemic(timevect,I0,npatches,onset_fixed,onset_thr,flag1,rs1,ps1,as1,Ks1,odeOptions)
%SIMULATESUBEPIDEMIC Shared event-driven simulation for fit and forecast paths.
%   Initial states and incidence conventions are unchanged from the initial-
%   condition patch. Every integration segment uses a fixed, explicit active
%   mask. At the first upward threshold crossing, activate the corresponding
%   patch(es) OUTSIDE the derivative and restart from the event state.
%
%   Returns exactly TIMEVECT(:), including for a two-point grid. Event times
%   are not inserted into the incidence grid. States remain continuous across
%   events; no seeds are added or removed at an activation. The full requested
%   grid must start at the calibration initial time (not a forecast-only tail).
%
%   Optional odeOptions retains ode15s tolerances and step controls. Events is
%   reserved for activation; nonempty custom Events, Mass, Jacobian, and
%   Vectorized='on' are rejected rather than silently misapplied to a changing
%   activation mask. The default solver remains ode15s.
%
%   Optional fifth output ACTIVATION contains:
%     active           final logical activation mask;
%     time             initial time for initially active patches, crossing
%                      time for triggered patches, NaN for never activated;
%     predecessorCount predecessor state at triggering, otherwise NaN;
%     initialActive    mask BEFORE threshold-at-start processing;
%     nSegments        number of smooth-segment ode15s calls.
%
%   Compatibility: legacy globals invasions, timeinvasions, Cinvasions are
%   output-only summaries, published AFTER a successful solve. They are never
%   read by the derivative, events, or integration control. Initially active
%   and never-active entries retain legacy zero times/counts in those globals;
%   use ACTIVATION to distinguish them. On entry globals are cleared, so a
%   failed call cannot leave a previous successful activation summary behind.
%   The surrounding optimizer still has other globals; this is not a claim
%   that the whole toolbox is thread-safe.
%
%   See also initializeSubepidemicState, modifiedLogisticGrowthPatch,
%            subepidemicActivationEvents, subepidemicIncidence

global invasions timeinvasions Cinvasions
invasions=[]; timeinvasions=[]; Cinvasions=[];

if nargin<11 || isempty(odeOptions)
    odeOptions=odeset();
end
[IC,active0,legacyTimes,legacyCounts]= ...
    initializeSubepidemicState(I0,npatches,onset_fixed);
if ~(isnumeric(timevect) && isreal(timevect) && isvector(timevect) && ...
        ~isempty(timevect) && all(isfinite(timevect(:))) && ...
        all(diff(timevect(:))>0))
    error('SubEpiPredict:InvalidTimeGrid', ...
        'timevect must be a finite, real, strictly increasing time vector.');
end
if ~(isnumeric(onset_thr) && isreal(onset_thr) && isscalar(onset_thr) && ...
        isfinite(onset_thr) && onset_thr>=0)
    error('SubEpiPredict:InvalidThreshold','onset_thr must be finite and nonnegative.');
end
if ~(isnumeric(flag1) && isreal(flag1) && isscalar(flag1) && ...
        isfinite(flag1) && flag1==fix(flag1) && flag1>=0 && flag1<=5)
    error('SubEpiPredict:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
end
if any([numel(rs1),numel(ps1),numel(as1),numel(Ks1)]~=npatches)
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Each growth-parameter vector must contain exactly npatches entries.');
end
parameters={rs1,ps1,as1,Ks1};
for j=1:numel(parameters)
    v=parameters{j};
    if ~(isnumeric(v) && isreal(v) && isvector(v) && all(isfinite(v(:))))
        error('SubEpiPredict:InvalidGrowthParameters', ...
            'Growth parameters must be finite, real numeric vectors.');
    end
end
if ~(isstruct(odeOptions) && isscalar(odeOptions))
    error('SubEpiPredict:InvalidOdeOptions','odeOptions must be an odeset structure.');
end
if ~isempty(odeget(odeOptions,'Events',[]))
    error('SubEpiPredict:ReservedEventsOption', ...
        'The simulator manages Events internally for subepidemic activation.');
end
if ~isempty(odeget(odeOptions,'Mass',[])) || ...
        ~isempty(odeget(odeOptions,'Jacobian',[])) || ...
        strcmpi(odeget(odeOptions,'Vectorized','off'),'on')
    error('SubEpiPredict:UnsupportedOdeOption', ...
        'Custom Mass/Jacobian and Vectorized=''on'' are not supported by this simulator.');
end

t=double(timevect(:));
x=nan(numel(t),npatches);
x(1,:)=IC.';
active=logical(active0);
initialActive=active;
activationTimes=nan(npatches,1);
activationTimes(active)=t(1);
activationCounts=nan(npatches,1);
currentTime=t(1);
currentState=double(IC);
nSegments=0;

% Preserve the literal original rule: x(j-1)>=threshold, including when a
% predecessor is a dormant one-case seed. No new predecessor-active condition
% is imposed. In particular, a zero threshold activates all eligible patches
% at t(1), without redistributing the asynchronous initial states.
activateAtBoundary();

while currentTime<t(end)
    nSegments=nSegments+1;
    if nSegments>npatches
        error('SubEpiPredict:ActivationProgressFailure', ...
            'Too many segments: every restart must activate a new patch.');
    end
    segmentActive=active; % anonymous callbacks capture this fixed VALUE
    watched=find(~segmentActive);
    watched=watched(watched>=2);
    rhs=@(tt,xx)modifiedLogisticGrowthPatch(tt,xx,rs1,ps1,as1,Ks1, ...
        npatches,onset_thr,flag1,segmentActive);
    if isempty(watched)
        segmentOptions=odeset(odeOptions,'Events',[]);
    else
        eventFcn=@(tt,xx)subepidemicActivationEvents(tt,xx,onset_thr,watched);
        segmentOptions=odeset(odeOptions,'Events',eventFcn);
    end
    segmentGrid=[currentTime;t(t>currentTime)];
    sol=ode15s(rhs,segmentGrid,currentState,segmentOptions);

    hasEvent=isfield(sol,'xe') && ~isempty(sol.xe);
    triggered=[];
    if hasEvent
        % MATLAB can report a first-step terminal event as nonterminal.
        % Always cut at the EARLIEST recorded event, discard any later old-
        % mask continuation, and restart from its event state, not sol.y(end).
        eventTimes=sol.xe(:);
        [stopTime,firstEvent]=min(eventTimes);
        stopState=sol.ye(:,firstEvent);
        eventRoundoff=8*eps(max(1,abs(stopTime)));
        simultaneous=abs(eventTimes-stopTime)<=eventRoundoff;
        eventIndices=sol.ie(simultaneous);
        triggered=unique(watched(eventIndices));
    else
        stopTime=sol.x(end);
        stopState=sol.y(:,end);
        if stopTime~=t(end)
            error('SubEpiPredict:IncompleteIntegration', ...
                'ode15s stopped before the requested endpoint without an activation event.');
        end
    end
    if ~(isfinite(stopTime) && stopTime>=currentTime && stopTime<=t(end)) || ...
            ~isreal(stopState) || any(~isfinite(stopState(:)))
        error('SubEpiPredict:InvalidSegment','Invalid event/segment endpoint.');
    end

    % Dense output is evaluated only within this smooth segment. Never
    % interpolate across a switch, insert event rows, or advance by an epsilon.
    rows=find(t>currentTime & t<=stopTime);
    if ~isempty(rows)
        x(rows,:)=deval(sol,t(rows).').';
    end
    currentTime=stopTime;
    currentState=stopState(:);

    if hasEvent
        oldCount=sum(active);
        activatePatches(triggered);
        activateAtBoundary(); % simultaneous or already-satisfied conditions
        if sum(active)<=oldCount
            error('SubEpiPredict:ActivationProgressFailure', ...
                'An event did not activate a previously dormant patch.');
        end
    else
        % A crossed threshold without a recorded event must not be silently
        % treated as a late activation at the final observation.
        dormant=find(~active);
        dormant=dormant(dormant>=2);
        if any(currentState(dormant-1)>=onset_thr)
            error('SubEpiPredict:MissedActivationEvent', ...
                'A threshold was passed without a located activation event.');
        end
    end
end

if ~isreal(x) || any(~isfinite(x(:)))
    error('SubEpiPredict:InvalidTrajectory', ...
        'Simulation must return a finite, real state at every requested time.');
end
if nargout>2
    [totinc,patchIncidence]=subepidemicIncidence(x,onset_fixed);
end
if nargout>4
    activation=struct('active',active,'time',activationTimes, ...
        'predecessorCount',activationCounts,'initialActive',initialActive, ...
        'nSegments',nSegments);
end
% Preserve existing sum(invasions) callers without using globals as inputs.
invasions=double(active);
timeinvasions=legacyTimes;
Cinvasions=legacyCounts;

    function activateAtBoundary()
        dormant=find(~active);
        dormant=dormant(dormant>=2);
        eligible=dormant(currentState(dormant-1)>=onset_thr);
        activatePatches(eligible);
    end

    function activatePatches(indices)
        indices=indices(:);
        indices=indices(~active(indices));
        if isempty(indices)
            return
        end
        active(indices)=true;
        activationTimes(indices)=currentTime;
        activationCounts(indices)=currentState(indices-1);
        legacyTimes(indices)=currentTime;
        legacyCounts(indices)=currentState(indices-1);
    end
end
