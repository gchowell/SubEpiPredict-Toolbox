function [P,fval,fitInfo] = fitSubepidemicModel(context,z0,LB,UB,isLog,starts,numstartpoints,settings)
%FITSUBEPIDEMICMODEL MultiStart fit of one candidate in a scaled search space.
%   [P,fval,fitInfo] = fitSubepidemicModel(context,z0,LB,UB,isLog,starts,numstartpoints,settings)
%
%   Port of fit_model.m of the sub-epidemic SEIR toolbox (revision 3.1) to
%   the n-sub-epidemic model. It replaces the MultiStart blocks of
%   fittingModifiedLogisticFunctionPatchABC.m and ...Multiple.m.
%
%   INPUT
%     context         fit context of subepidemicObjective.m (data, times,
%                     initial value, number of sub-epidemics, onset rule,
%                     threshold, kernel, method, ODE options)
%     z0, LB, UB      seed and bounds, layout [r, p, a, K, alpha, d] with
%                     4*n+2 entries (subepidemicSearchSpace.m). A parameter
%                     with LB == UB is fixed at that value.
%     isLog           logical, same length: coordinates searched on a log scale
%     starts          [] or a matrix of start points, one per row in original
%                     units, used for the FIRST MultiStart pass (the screen
%                     builds one set per number of sub-epidemics and reuses
%                     it for every threshold). With [] the function draws
%                     numstartpoints Latin-hypercube points, jittered copies
%                     of z0 and z0 itself.
%     numstartpoints  Latin-hypercube points per generated start set
%     settings        subepidemicFitSettings.m (optional)
%
%   OUTPUT
%     P       fitted vector in original units (NaN when the fit failed)
%     fval    objective at P in the legacy units (SSE, negative
%             log-likelihood without its data constant, or SAD): the number
%             getAICc expects. Inf when the fit failed.
%     fitInfo diagnostics, see below. fitInfo.success must be checked.
%
%   WHAT THE OPTIMIZER SEES
%   1. Only the free parameters, each mapped to [0, 1] (subepidemicUnitMap).
%      Before, fmincon received the whole vector, fixed entries included,
%      with K up to 1e10 beside p in [0, 1]. Its relative step test was then
%      decided by K, and local runs stopped after one or two iterations with
%      a positive exit flag far from any minimum. In the bundled COVID-19
%      negative-binomial run K stayed within 83 of its seed 84,101 and the
%      four top-ranked fits were 18 to 32 log-likelihood units above the
%      minimum next to them.
%   2. The objective in log-likelihood units, divided by the number of
%      observations so it is of order 1:
%          least squares        (n/2) * log(SSE)
%          absolute deviations   n * log(SAD)
%          Poisson, neg. bin.    legacy objective + sum(gammaln(floor(y)+1))
%      These are monotone transformations of the legacy objective (the last
%      one adds a constant), so the minimizer is unchanged. The scale no
%      longer depends on how good the seed is, and tolerances and
%      diagnostics have the same meaning for every estimation method.
%      (fit_model.m divides by the objective at the seed instead.)
%   An invalid parameter vector (failed integration, zero trajectory) gets a
%   penalty of 1e100; every local minimum is evaluated again in original
%   units and must be valid before it is accepted.
%
%   DIAGNOSTICS in fitInfo (a positive exit flag is not a minimum)
%     success, message, exitflag, attempts
%     nActive               sub-epidemics active by the end of the calibration
%                           grid at P. The screen keeps a candidate only when
%                           this equals the number of sub-epidemics fitted.
%     startAgreement        fraction of converged local runs whose objective is
%                           within settings.agreementTol log-likelihood units
%                           of the best converged run. Near 1: the minimum was
%                           found from almost every start. Small: the surface
%                           is multimodal and the number of starts matters.
%     nLocalRuns, nConvergedRuns
%     gradient              projected central-difference gradient at P, in
%                           log-likelihood units per unit coordinate (one unit
%                           = the full range of that parameter); components
%                           pushing into an active bound are zero
%     gradientInfNorm       its largest absolute component. In the tests of
%                           this patch a fit at its minimum scored below 1 and
%                           the stalled fits of the previous version about
%                           2,000; the polish gate is 10.
%     polish                struct: enabled, run, improved, exitflag,
%                           confirmExitflag, objectiveBefore, objectiveAfter,
%                           gradientBefore, gradientTol, budget. The
%                           Nelder-Mead polish runs only when gradientInfNorm
%                           exceeds settings.polishGradientTol.
%     atLowerBound, atUpperBound   indices (into P) of free parameters at a bound
%     fallback              when success is false: struct P, fval, nActive of
%                           the best VALID point any local run ended at, or [].
%                           The bootstrap uses it to keep its legacy policy of
%                           retaining unresolved refits; the screen does not.
%     freeIndex, nFree, isLog, optimizerOptions, useParallel, objectiveUnits
%
%   Random numbers come from the global stream (start points).
%
%   See also subepidemicObjective, subepidemicSearchSpace, subepidemicUnitMap,
%            subepidemicStartPoints, subepidemicFitSettings

if nargin<6, starts = []; end
if nargin<7 || isempty(numstartpoints), numstartpoints = 0; end
if nargin<8 || isempty(settings), settings = subepidemicFitSettings(); end
validateattributes(numstartpoints,{'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative'});
required = {'timevect','ydata','I0','npatches','onset_fixed','onset_thr','flag1','method1','odeOptions'};
for k = 1:numel(required)
    if ~isfield(context,required{k})
        error('SubEpiPredict:FitContext','The fit context has no field ''%s''.',required{k});
    end
end
if ~any(context.method1==[0 1 3 4 5 6])
    error('SubEpiPredict:InvalidMethod','method1 must be 0, 1, 3, 4, 5 or 6.');
end
nz = 4*full(double(context.npatches))+2;
if any([numel(z0),numel(LB),numel(UB),numel(isLog)]~=nz)
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Seed, bounds and scale flags must all have 4*npatches+2 = %d entries.',nz);
end
map = subepidemicUnitMap('make',struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog));
z0 = map.z0;
nFree = map.nFree;
ydata = double(context.ydata(:));
nObs = numel(ydata);
if nObs<2 || any(~isfinite(ydata)) || ~isreal(ydata)
    error('SubEpiPredict:InvalidData','The fitted series must hold at least two finite real values.');
end

% Objective in log-likelihood units. The data constant restores the term the
% count objectives drop, so their value is a full negative log-likelihood.
dataConstant = 0;
if any(context.method1==[1 3 4 5])
    dataConstant = sum(gammaln(floor(ydata)+1));
end
natU = @(u) localNatural(subepidemicUnitMap('toZ',u,map),context,nObs,dataConstant);
fs = @(u) natU(u)/nObs;                                                % what fmincon minimizes

P = NaN(1,nz); fval = Inf;
fitInfo = struct('success',false,'message','','exitflag',NaN,'attempts',0, ...
    'attemptExitflags',[],'nActive',0, ...
    'freeIndex',map.free,'nFree',nFree,'isLog',map.isLog, ...
    'objectiveUnits','log-likelihood units; fmincon minimizes them divided by the number of observations', ...
    'useParallel',false,'nLocalRuns',0,'nConvergedRuns',0,'startAgreement',NaN, ...
    'gradient',[],'gradientInfNorm',NaN, ...
    'atLowerBound',[],'atUpperBound',[],'fallback',[]);
fitInfo.attemptMessages = {};
fitInfo.optimizerOptions = struct();
fitInfo.polish = struct('enabled',settings.polish,'run',false,'improved',false, ...
    'exitflag',NaN,'confirmExitflag',NaN,'objectiveBefore',NaN,'objectiveAfter',NaN, ...
    'gradientBefore',NaN,'gradientTol',settings.polishGradientTol,'budget',settings.polishBudget);

accepted = false; chosenDetail = [];
fallback = [];

if nFree==0
    % nothing to estimate: evaluate the fixed vector
    [value,detail] = subepidemicObjective(z0,context);
    if detail.success
        P = z0; fval = value; chosenDetail = detail; accepted = true;
        fitInfo.exitflag = 1;
    else
        fitInfo.message = detail.message;
    end
else
    % ---- optimizer settings (fit_model.m defaults, then the user's overrides) --
    options = optimoptions('fmincon','Algorithm','sqp', ...
        'StepTolerance',1e-10,'FunctionTolerance',1e-10, ...
        'OptimalityTolerance',1e-8,'FiniteDifferenceType','forward', ...
        'FiniteDifferenceStepSize',1e-6,'MaxFunctionEvaluations',20000, ...
        'MaxIterations',20000,'Display','off');
    keys = fieldnames(settings.optimOptions);
    for j = 1:numel(keys)
        options = optimoptions(options,keys{j},settings.optimOptions.(keys{j}));
    end
    fdStep = options.FiniteDifferenceStepSize;
    if ~isscalar(fdStep), fdStep = fdStep(1); end
    fitInfo.optimizerOptions = struct('Algorithm',options.Algorithm, ...
        'StepTolerance',options.StepTolerance,'FunctionTolerance',options.FunctionTolerance, ...
        'OptimalityTolerance',options.OptimalityTolerance, ...
        'FiniteDifferenceStepSize',fdStep,'FiniteDifferenceType',options.FiniteDifferenceType, ...
        'MaxFunctionEvaluations',options.MaxFunctionEvaluations,'MaxIterations',options.MaxIterations);
    problem = createOptimProblem('fmincon','objective',fs, ...
        'x0',subepidemicUnitMap('toU',z0,map), ...
        'lb',zeros(1,nFree),'ub',ones(1,nFree),'options',options);

    % ---- MultiStart, serial or parallel ---------------------------------------
    % The objective captures its context in the handle, so workers need no
    % globals. (simulateSubepidemic still publishes its legacy output-only
    % globals; nothing reads them during a fit.)
    if settings.useParallel
        try
            if isempty(gcp('nocreate')), parpool; end
            ms = MultiStart('Display','off','UseParallel',true);
            fitInfo.useParallel = true;
        catch ME
            subepidemicWarning('SubEpiPredict:ParallelUnavailable', ...
                'Parallel pool unavailable (%s); running MultiStart serially.',ME.message);
            ms = MultiStart('Display','off','UseParallel',false);
        end
    else
        ms = MultiStart('Display','off','UseParallel',false);
    end

    nJitter = min(3,max(1,numstartpoints));
    for attempt = 1:settings.fitMaxAttempts
        if attempt==1 && ~isempty(starts)
            startSet = double(full(starts));
            if size(startSet,2)~=nz || ~isreal(startSet) || any(~isfinite(startSet(:)))
                error('SubEpiPredict:ParameterDimensionMismatch', ...
                    'Start points must be finite, one per row, with %d columns.',nz);
            end
        else
            startSet = subepidemicStartPoints(map,numstartpoints,nJitter,z0);
        end
        U0 = unique(subepidemicUnitMap('toU',startSet,map),'rows','stable');
        fitInfo.attempts = attempt;
        try
            [candidate,value,flag,output,minima] = run(ms,problem,CustomStartPointSet(U0));
        catch ME
            flag = -10; candidate = []; value = Inf; minima = [];
            output = struct('message',sprintf('%s: %s',ME.identifier,ME.message));
        end
        fitInfo.attemptExitflags(end+1) = flag;
        fitInfo.attemptMessages{end+1} = output.message;
        fitInfo.exitflag = flag;
        fitInfo.nLocalRuns = size(U0,1);
        fitInfo.nConvergedRuns = 0;

        % Start agreement: how many converged local runs reproduce the best one.
        if ~isempty(minima)
            fvs = [minima.Fval]*nObs;                      % log-likelihood units
            okRun = ([minima.Exitflag]>0) & isfinite(fvs);
            fitInfo.nConvergedRuns = sum(okRun);
            if any(okRun)
                bestOK = min(fvs(okRun));
                fitInfo.startAgreement = sum(okRun & fvs<=bestOK+settings.agreementTol)/sum(okRun);
            end
        end

        % A penalty plateau can appear to converge. Evaluate each local minimum
        % again in original units and accept the best one that is valid.
        candidates = {}; values = []; flags = [];
        if flag>0 && ~isempty(candidate) && isfinite(value)
            candidates{end+1} = subepidemicUnitMap('toZ',candidate,map); %#ok<AGROW>
            values(end+1) = value; flags(end+1) = flag; %#ok<AGROW>
        end
        for j = 1:numel(minima)
            if minima(j).Exitflag>0 && isfinite(minima(j).Fval)
                candidates{end+1} = subepidemicUnitMap('toZ',minima(j).X,map); %#ok<AGROW>
                values(end+1) = minima(j).Fval; %#ok<AGROW>
                flags(end+1) = minima(j).Exitflag; %#ok<AGROW>
            end
        end
        [~,order] = sort(values);
        for j = order
            trial = candidates{j}(:).';
            [value,detail] = subepidemicObjective(trial,context);
            if detail.success
                P = trial; fval = value; chosenDetail = detail;
                fitInfo.exitflag = flags(j); accepted = true;
                break
            end
            fitInfo.message = detail.message;
        end
        if accepted, break; end

        % Not accepted: remember the best valid end point for the fallback record.
        if ~isempty(candidate) && isreal(candidate) && all(isfinite(candidate(:))) && numel(candidate)==nFree
            trial = subepidemicUnitMap('toZ',candidate,map);
            [value,detail] = subepidemicObjective(trial,context);
            if detail.success && (isempty(fallback) || value<fallback.fval)
                fallback = struct('P',trial,'fval',value,'nActive',detail.nActive);
            end
        end
        % An unsuccessful pass is retried with fresh start points.
    end

    % ---- stationarity check, then (only if needed) the polish ----------------
    if accepted
        [g,gInf] = localGradient(natU,subepidemicUnitMap('toU',P,map),fdStep);
        fitInfo.gradient = g; fitInfo.gradientInfNorm = gInf;
        fitInfo.polish.gradientBefore = gInf;
    end
    if accepted && settings.polish && ...
            ~(isfinite(fitInfo.gradientInfNorm) && fitInfo.gradientInfNorm<=settings.polishGradientTol)
        fitInfo.polish.run = true; fitInfo.polish.objectiveBefore = fval;
        uP = subepidemicUnitMap('toU',P,map);
        budget = ceil(settings.polishBudget*nFree);
        optsNM = optimset('TolX',1e-6,'TolFun',1e-6,'Display','off', ...
            'MaxFunEvals',budget,'MaxIter',budget);
        try
            [uRaw,~,polishFlag] = fminsearch(@(u) localBoxPenalty(u,fs),uP,optsNM);
            fitInfo.polish.exitflag = polishFlag;
            uTry = min(max(uRaw(:).',0),1);
            zTry = subepidemicUnitMap('toZ',uTry,map);
            [vTry,dTry] = subepidemicObjective(zTry,context);   % validity check in original units
            if dTry.success && vTry<fval
                P = zTry; fval = vTry; chosenDetail = dTry; fitInfo.polish.improved = true;
                % A polish that improved but did not converge (budget used up, or
                % an optimum on a bound) is handed to a bounded sqp run from the
                % polished point; its exit flag describes the returned vector.
                if polishFlag~=1 || any(uRaw(:).'<0) || any(uRaw(:).'>1)
                    try
                        [uC,vC,flagC] = fmincon(fs,uTry,[],[],[],[], ...
                            zeros(1,nFree),ones(1,nFree),[],options); %#ok<ASGLU>
                        fitInfo.polish.confirmExitflag = flagC;
                        zC = subepidemicUnitMap('toZ',uC,map);
                        [vC2,dC] = subepidemicObjective(zC,context);
                        if flagC>0 && dC.success && vC2<=fval
                            P = zC; fval = vC2; chosenDetail = dC; fitInfo.exitflag = flagC;
                        end
                    catch
                        % keep the polished point
                    end
                end
            end
        catch
            % keep the sqp point if the polish fails for any reason
        end
        fitInfo.polish.objectiveAfter = fval;
        if fitInfo.polish.improved        % check stationarity again at the returned point
            [g,gInf] = localGradient(natU,subepidemicUnitMap('toU',P,map),fdStep);
            fitInfo.gradient = g; fitInfo.gradientInfNorm = gInf;
        end
    end
end

if accepted
    fitInfo.success = true;
    fitInfo.nActive = chosenDetail.nActive;
    fitInfo.message = 'Accepted a converged local minimum with a valid trajectory.';
    if nFree>0
        uP = subepidemicUnitMap('toU',P,map);
        fitInfo.atLowerBound = map.free(uP<=1e-6);
        fitInfo.atUpperBound = map.free(uP>=1-1e-6);
    end
else
    % The seed is the last resort of the fallback record, never an accepted fit.
    [value,detail] = subepidemicObjective(z0,context);
    if detail.success && (isempty(fallback) || value<fallback.fval)
        fallback = struct('P',z0,'fval',value,'nActive',detail.nActive);
    end
    fitInfo.fallback = fallback;
    if isempty(fitInfo.message)
        fitInfo.message = ['No converged, valid local minimum. Seed diagnostic: ',detail.message];
    end
    P(:) = NaN; fval = Inf;
end
end

% <============================================================================>
% < Local functions                                                            >
% <============================================================================>
function value = localNatural(z,context,nObs,dataConstant)
% Objective in log-likelihood units; 1e100 for an invalid parameter vector.
[v,detail] = subepidemicObjective(z,context);
if ~detail.success
    value = 1e100;
    return
end
switch context.method1
    case 0
        value = 0.5*nObs*log(max(v,realmin));
    case 6
        value = nObs*log(max(v,realmin));
    otherwise
        value = v+dataConstant;
end
if ~isfinite(value)
    value = 1e100;
end
end

function [g,gInf] = localGradient(natU,uP,h)
% Projected central-difference gradient in unit coordinates, log-likelihood
% units. One-sided at a bound. Components that push into an active bound are
% set to zero. Central differences have no second-order bias, so the check
% does not fail at a converged fit because of the step length.
g = NaN(1,numel(uP)); gInf = NaN;
try
    uP = uP(:).'; g = zeros(1,numel(uP)); fP = [];
    for j = 1:numel(uP)
        e = zeros(1,numel(uP)); e(j) = h;
        if uP(j)+h<=1 && uP(j)-h>=0
            g(j) = (natU(uP+e)-natU(uP-e))/(2*h);
        else
            if isempty(fP), fP = natU(uP); end
            if uP(j)+h<=1
                g(j) = (natU(uP+e)-fP)/h;
            else
                g(j) = (fP-natU(uP-e))/h;
            end
        end
    end
    g(uP<=1e-8 & g>0) = 0;          % at the lower bound, descent would leave the box
    g(uP>=1-1e-8 & g<0) = 0;        % at the upper bound, likewise
    gInf = norm(g,inf);
catch
    % diagnostics only
end
end

function v = localBoxPenalty(u,fs)
% Nelder-Mead objective: the scaled objective inside the unit box, a penalty outside.
u = u(:).';
if any(u<0) || any(u>1)
    v = 1e100; return
end
v = fs(u);
end
