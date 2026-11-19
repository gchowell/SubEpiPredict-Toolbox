% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [RMSES,PS,npatches,onset_thr,P]=fittingModifiedLogisticFunctionPatchABC(datafilename1,data1,DT,epidemic_period,M,flagX,numstartpoints)

global flag1 method1 timevect ydata

global I0 npatches onset_thr

flag1=flagX;

close all

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global dist1
global factor1

global smoothfactor1
global calibrationperiod1

global LBe UBe

% <=============================================================================================>
% <========= Set bounds for parameters associated with error structure (alpha and d) ===========>
% <=============================================================================================>

switch method1
    
    case 0
        LBe=[0 0];
        UBe=[0 0];
    case 1
        LBe=[0 0];
        UBe=[0 0];
    case 2
        LBe=[0 0];
        UBe=[0 0];
    case 3
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 4
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 5
        LBe=[10^-8 0.2];
        UBe=[10^5 10^2];
end

% <=============================================================================================>
% <========= Fitting layer (revision 2026-11: scaled search space) =============================>
% <=============================================================================================>
% Every candidate is now fitted by fitSubepidemicModel.m: only the free
% parameters are optimized, each mapped to [0,1] (r, K and alpha on a log
% scale), the objective is minimized in log-likelihood units, every local
% minimum is evaluated again before it is accepted, and each fit returns
% diagnostics (start agreement, gradient, polish). Bounds and seeds are the
% legacy ones unless options_fit_settings.m changes them
% (subepidemicSearchSpace.m). The ODE tolerance is set in
% subepidemicOdeOptions.m. RMSES, PS and thresholdDiagnostics keep their
% layout; fitDiagnostics and candidateFitInfo are new and are saved with them.
% The globals LBe and UBe above are kept for old scripts and are no longer
% read by the fit.

fitSettings=subepidemicFitSettings();
odeOptions=subepidemicOdeOptions();

% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>

%data1=load(strcat('./output/',datafilename1));

data1=data1(epidemic_period,:);

data1(:,2)=data1(:,2);

I0=data1(1,2); % initial condition

if I0==0
    data1=data1(2:end,:);
end

data=data1(:,2);

% <==============================================================================>
% <============ Set time vector (timevect) and initial condition (I0) ===========>
% <==============================================================================>

timevect=(data1(:,1));

I0=data(1); % initial condition

% <==============================================================================>
% <===================== Set initial parameter guesses ==========================>
% <==============================================================================>

rs1=zeros(1,npatches_fixed)+0.1;
ps1=zeros(1,npatches_fixed)+0.9;
Ks1=zeros(1,npatches_fixed)+sum(data1(:,2));

as1=ones(1,npatches_fixed);

for j=1:npatches_fixed
    
    if flag1==3 | flag1==4 | flag1==5 % Logistic model or Richards model (p=1)
        ps1(j)=1;
    else
        ps1(j)=0.9;
    end
    
    if flag1==5
        
        rs1(j)=1-I0/Ks1(j);
        
        as1(j)=rs1(j)/log(Ks1(j)/I0);
        
    end
    
end

% <==============================================================================>
% <================= Set range of C_thr values (onset_thrs) =====================>
% <==============================================================================>

%cum1=sum(smooth(data1(:,2),smoothfactor1));
%onset_thrs=unique(cumsum(smooth(data1(:,2),smoothfactor1)));
%index2=find(onset_thrs<=0.99*cum1);
%onset_thrs=onset_thrs(index2)';

% equal-witdh discretization of C_thr
cumcurve1=cumsum(smoothdata(data1(:,2),'movmean',smoothfactor1));

onset_thrs=linspace(cumcurve1(1),cumcurve1(end),length(data1(:,2)));

onset_thrs=[0 onset_thrs(1:end-1)];

% <=============================================================================>
% <===== Set range of the possible number of subepidemics (1:npatches_fixed)====>
% <=============================================================================>

npatchess=1:1:npatches_fixed;

if onset_fixed==1 | (length(npatchess)==1 & npatchess(1)==1)
    
    onset_thrs=0;
    
end

onset_thrs2=onset_thrs;

RMSES=sparse(1000,3);

PS=sparse(1000,npatches_fixed*4+2);

count1=1;

% Diagnostics are in search order, including skipped partial-activation fits.
% Keep RMSES (3 columns) and PS (4*npatches_fixed+2 columns) unchanged.
thresholdDiagnosticsColumns={'nCandidate','onset_thr','nActive','retained','fval'};
thresholdDiagnostics=nan(numel(npatchess)*numel(onset_thrs2),5);
nEvaluated=0;
candidatesTooLarge=zeros(0,3); % [sub-epidemics, parameters, data points] of sizes left out (revision h)

% One row per fitted candidate, in the same search order as thresholdDiagnostics.
% candidateFitInfo{k} is the full fitInfo structure of row k.
fitDiagnosticsColumns={'nCandidate','onset_thr','success','exitflag','startAgreement', ...
    'gradientInfNorm','polishRun','polishImproved','nLocalRuns','nConvergedRuns','nAtBound'};
fitDiagnostics=nan(numel(npatchess)*numel(onset_thrs2),numel(fitDiagnosticsColumns));
candidateFitInfo=cell(numel(npatchess)*numel(onset_thrs2),1);

% <====================================================================================>
% <==== Evaluate AICc across models with different number of subepidemics and C_thr ==>
% <====================================================================================>

% Series the objective is fitted to (revision 2026-11-16). Count likelihoods
% (methods 1, 3, 4, 5) are fitted to the raw counts; least squares and
% absolute deviations to the moving average, as before. smoothfactor1 still
% shapes the threshold grid above for every method. See subepidemicFitSeries.m.
[ydata,fitSeriesInfo]=subepidemicFitSeries(data,smoothfactor1,method1,fitSettings);

if fitSeriesInfo.smoothedCountLikelihood
    % Only reached when countLikelihoodOnRawData = false in options_fit_settings.m.
    subepidemicWarning('SubEpiPredict:SmoothedCountLikelihood', ...
        ['method1 = %d fits a count likelihood to a %g-point moving average. The fitted dispersion ', ...
         'will be too small and prediction intervals too narrow for the raw counts.'],method1,smoothfactor1);
elseif any(method1==[1 3 4 5]) && smoothfactor1>1
    fprintf(['method1 = %d is fitted to the raw counts. smoothfactor1 = %g is used only for the grid of ', ...
        'onset thresholds.\n'],method1,smoothfactor1);
end

for npatches2=[npatchess]
    
    % The optimizer dimension is fixed for this entire threshold sweep.
    % nActive is a separate observation, never a replacement for nCandidate.
    nCandidate=npatches2;
    npatches=nCandidate;

    % Revision h. AICc is defined only when the points fitted exceed the
    % number of parameters by more than one (see getAICc.m). A candidate that
    % is too large for this series is not fitted: its AICc would be infinite
    % or, worse, lowered by a negative small-sample term.
    kCandidate=get_nparams(method1,dist1,nCandidate,flag1,1,onset_fixed);
    if length(ydata)<=kCandidate+1
        candidatesTooLarge(end+1,:)=[nCandidate kCandidate length(ydata)]; %#ok<AGROW>
        subepidemicWarning('SubEpiPredict:CandidateTooLarge', ...
            ['Models with %d sub-epidemic(s) have %d parameters, too many for the %d data points ', ...
             '(AICc needs more than %d). They are left out of the screen.'], ...
            nCandidate,kCandidate,length(ydata),kCandidate+1);
        continue
    end
    
    if (onset_fixed==1 | npatches==1)
        onset_thrs=0;
    else
        onset_thrs=onset_thrs2;
        % A threshold that the preceding sub-epidemic already meets at time
        % zero starts every sub-epidemic at once, exactly as threshold 0
        % does: the first sub-epidemic starts at I0 and each later one at 1.
        % The grid's second value, cumcurve1(1), is usually such a threshold
        % (always when smoothfactor1 = 1). Fitting it would rank the same
        % model twice and put it twice into the top-ranked ensemble.
        if nCandidate==2
            levelAtStart=I0;
        else
            levelAtStart=min(I0,1);
        end
        metAtStart=find(onset_thrs<=levelAtStart);
        if numel(metAtStart)>1
            % All of these are one model; keep the first (threshold 0).
            sameAsZero=metAtStart(2:end);
            fprintf('Onset thresholds dropped as duplicates of threshold %g (already met at time zero): %s\n', ...
                onset_thrs(metAtStart(1)),mat2str(onset_thrs(sameAsZero),6));
            onset_thrs(sameAsZero)=[];
        end
    end
    
    % <================================================================================================>
    % <=========================== Search space, seed and start points ===============================>
    % <================================================================================================>

    % Legacy bounds and seed unless options_fit_settings.m overrides them.
    % isLog marks the coordinates searched on a log scale (r, K, alpha).
    [LB,UB,isLog,z0]=subepidemicSearchSpace(nCandidate,data,flag1,method1,fitSettings);
    expectedDimension=4*nCandidate+2;
    if any([numel(z0),numel(LB),numel(UB)]~=expectedDimension)
        error('SubEpiPredict:ParameterDimensionMismatch', ...
            'Seed and bounds must have 4*nCandidate+2 entries.');
    end
    z=z0; % legacy name of the seed, kept in the saved workspace

    % Start points ONCE per nCandidate, reused for every onset threshold:
    % numstartpoints Latin-hypercube points in unit coordinates (log-uniform
    % for r, K and alpha) and the seed. Rows are in original units.
    unitMap=subepidemicUnitMap('make',struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog));
    nStarts=max(0,floor(numstartpoints));
    starts_base=subepidemicStartPoints(unitMap,nStarts,0,z0);

    % Accepted solution of the previous threshold for this nCandidate
    % (continuation start); empty until one candidate has been retained.
    previousP=[];

    nloops=length(onset_thrs);
    
    RMSES2=sparse(1000,3);
    count2=1;

    for onset_thr=onset_thrs
        % The candidate dimension is fixed for the whole threshold sweep.
        npatches=nCandidate;

        % Explicit fit context: the objective reads nothing from globals.
        context=struct('timevect',timevect,'ydata',ydata,'I0',I0, ...
            'npatches',nCandidate,'onset_fixed',onset_fixed,'onset_thr',onset_thr, ...
            'flag1',flag1,'method1',method1,'odeOptions',odeOptions);

        % The solution of the previous threshold is an excellent start for
        % this one: for a fixed nCandidate the best fit changes smoothly with
        % the threshold. It is one more start point, never a replacement.
        starts=starts_base;
        nStartsThisFit=nStarts;
        if fitSettings.continuationStart && ~isempty(previousP)
            starts=[starts_base;previousP];
        elseif isempty(previousP) && fitSettings.firstThresholdStartFactor>1 && nStarts>0
            % Revision m. A fit without a solution of a previous threshold to
            % start from (the first threshold of the sweep, and any later one
            % until a candidate has been retained) depends on the random
            % start points alone. With 20 of them the search of the first
            % threshold of the COVID-19 example ended in a poor local optimum
            % for about 1 seed in 10, and a model with an evidence ratio of 2
            % dropped out of the ranking. Such a fit gets
            % firstThresholdStartFactor times the start points (default 3).
            nExtra=ceil((fitSettings.firstThresholdStartFactor-1)*nStarts);
            extraStarts=subepidemicStartPoints(unitMap,nExtra,0,z0);
            starts=[extraStarts(1:end-1,:);starts_base]; % the seed is already the last row of starts_base
            nStartsThisFit=nStarts+nExtra;
        end

        [P,fval,fitInfo]=fitSubepidemicModel(context,z0,LB,UB,isLog,starts,nStartsThisFit,fitSettings);

        nEvaluated=nEvaluated+1;
        candidateFitInfo{nEvaluated}=fitInfo;
        fitDiagnostics(nEvaluated,:)=[nCandidate onset_thr double(fitInfo.success) fitInfo.exitflag ...
            fitInfo.startAgreement fitInfo.gradientInfNorm double(fitInfo.polish.run) ...
            double(fitInfo.polish.improved) fitInfo.nLocalRuns fitInfo.nConvergedRuns ...
            numel(fitInfo.atLowerBound)+numel(fitInfo.atUpperBound)];

        % A candidate with no converged, valid local minimum is recorded and
        % skipped. It is never replaced by its seed.
        if ~fitInfo.success
            thresholdDiagnostics(nEvaluated,:)=[nCandidate onset_thr NaN 0 NaN];
            subepidemicWarning('SubEpiPredict:CandidateRejected', ...
                'Rejected n=%g, onset_thr=%g: %s',nCandidate,onset_thr,fitInfo.message);
            continue
        end

        % Activation is read from the accepted fit's own simulation.
        nActive=fitInfo.nActive;
        thresholdDiagnostics(nEvaluated,:)=[nCandidate onset_thr nActive ...
            double(nActive==nCandidate) fval];

        % Conservative policy: do not relabel or rescore a partially active
        % nCandidate fit as a smaller model. Every smaller dimension is fitted
        % independently by the outer loop, with its own bounds and starts.
        % A deliberately reduced candidate would need a separate refit.
        if nActive~=nCandidate
            continue
        end

        % Revision m: sub-epidemics that start together are interchangeable;
        % give them a fixed order (subepidemicOrderLabels.m). The curve and
        % the AICc do not change.
        P=subepidemicOrderLabels(P,nCandidate,onset_fixed,onset_thr,I0);

        previousP=P;

        AICc=getAICc(method1,dist1,nCandidate,flag1,1,fval,length(ydata),onset_fixed);

        RMSES(count1,:)=[nCandidate onset_thr AICc];
        PS(count1,1:expectedDimension)=P;
        count1=count1+1;

    end %onset
    
end %npatches

thresholdDiagnostics=thresholdDiagnostics(1:nEvaluated,:);
fitDiagnostics=fitDiagnostics(1:nEvaluated,:);
candidateFitInfo=candidateFitInfo(1:nEvaluated);
if count1==1 && ~isempty(candidatesTooLarge) && nEvaluated==0
    error('SubEpiPredict:TooFewDataPoints', ...
        ['No model can be ranked: the series has %d data points and the smallest model has %d parameters. ', ...
         'AICc needs more data points than parameters plus one.'], ...
        length(ydata),min(candidatesTooLarge(:,2)));
end
if count1==1
    error('SubEpiPredict:NoFullyActiveCandidates', ...
        'No candidate activated all of its fitted subepidemics during calibration.');
end

%RMSES(1:count1,:)
%pause



% <=============================================================================================>
% <======================== Sort the results by AICc (lowest to highest) =======================>
% <=============================================================================================>

RMSES=RMSES(1:count1-1,:);

PS=PS(1:count1-1,:);

[RMSES,index1]=sortrows(RMSES,[3 1]);

PS=PS(index1,:);


[RMSE1, index1]=min(RMSES(:,3));

npatches=RMSES(index1,1);

onset_thr=RMSES(index1,2);

AICc_best=RMSES(index1,3);


 %-->If we have a series of AICc  (or wSSE) values from N different models sorted from lowest (best model) 
 %to highest (worst model), I am wondering if we could define a proper threshold criterion to drop models with associated AICc (or wSSE) value greater than some threshold criteria.

% -->Let AICmin denote the minimun AIC from several models. 
% The quantity exp((AICmin − AICi)/2) is interpreted as the relative likelihood of model i. 
% We can set an alpha (e.g., 0.05), drop the models with exp((AICmin − AICi)/2) smaller than alpha, 
% and combine other models with weighted average, 
% where the weight is proportional to exp((AICmin − AICi)/2). This is another way to assign weights,
% compared to 1/SSE. 


AICmin=RMSES(1,3);

relativelik_i=exp((AICmin-RMSES(:,3))/2);

% -> Drop models with alpha below 0.05
%index2=find(relativelik_i>0.05);
%RMSES=RMSES(index2,:);
%relativelik_i=relativelik_i(index2);

if 0

    subplot(1,2,1)
    line1=plot(RMSES(:,3),'ko-');

    set(line1,'LineWidth',2)

    xlabel('Model i')
    ylabel('AICc')
    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    subplot(1,2,2)


    line1=plot(relativelik_i,'ko-');

    set(line1,'LineWidth',2)

    xlabel('Model i')
    ylabel('Relative likelihood')
    set(gca,'FontSize', 16);
    set(gcf,'color','white')

end



% <=============================================================================================>
% <============================ Get the best fit resultls ======================================>
% <=============================================================================================>


npatches=full(npatches); onset_thr=full(onset_thr); AICc_best=full(AICc_best);

P=full(PS(index1,1:npatches*4+2));

[rs_hat,ps_hat,as_hat,Ks_hat,alpha_hat,d_hat]= ...
    unpackSubepidemicParameters(P,npatches);

if method1==3
    
    dist1=3; % VAR=mean+alpha*mean;
    
    factor1=alpha_hat;
    
elseif method1==4
    
    dist1=4; % VAR=mean+alpha*mean^2;
    
    factor1=alpha_hat;
    
elseif method1==5
    
    dist1=5; % VAR=mean+alpha*mean^2;
    
    factor1=alpha_hat;
    
end


% Shared initialization, activation reset, and incidence conversion.
[~,x,totinc,patchIncidence,bestActivation]=simulateSubepidemic(timevect,I0,npatches,onset_fixed, ...
    onset_thr,flag1,rs_hat,ps_hat,as_hat,Ks_hat,odeOptions);

% The stored dimension, parameter vector, and AICc must stay paired. Do not
% silently shrink the selected model after ranking if reconstruction differs.
nActiveBest=sum(bestActivation.active);
if nActiveBest~=npatches
    error('SubEpiPredict:BestModelActivationMismatch', ...
        ['The selected %d-subepidemic model reconstructed with %d active patches. ' ...
         'The saved model must be refitted, not truncated after ranking.'], ...
        npatches,nActiveBest);
end


% <=============================================================================================>
% <================================= Plot best model fit ======================================>
% <=============================================================================================>

figure(100)
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
nexttile(1)
for j=1:npatches
    
    incidence1=patchIncidence(:,j);
    
    plot(timevect,incidence1)
    hold on
    
end



bestfit=totinc;

plot(timevect,totinc,'r')

hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('Incidence')

title('best fit')
fprintf('Best model of the screen: %d sub-epidemic(s), onset threshold %g.\n',full(npatches),full(onset_thr));


if method1==0
    
    if dist1==2  % calculate the overdispersion factor
        
        %     [coef,ns]=getMeanVarLinear(data,totinc,6);
        %
        %     if coef>0
        %         factor1=coef;
        %     else
        
        
        % estimate dispersion in data
        binsize1=7; %4
        
        [ratios,~]=getMeanVarianceRatio(data,binsize1,2);  % **
        
        %[ratios,~]=getMeanVarianceRatio(data,binsize1,1);
        
        index1=find(ratios(:,1)>0);
        
        factor1=mean(ratios(index1,1));
        
        factor1;
        
    end
    
end


% <=============================================================================================>
% <===================================  Save the results  ======================================>
% <=============================================================================================>

save(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

