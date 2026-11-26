% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [Phatss,npatches,onset_thr,curves,bestfit,data1,P0,AICc_best,factor1,d,bootstrapDiagnostics]=fittingModifiedLogisticFunctionMultiple(RMSES,PS,data1,DT,epidemic_period,M,flagX,numstartpoints,rank1)

global flag1 timevect ydata yfit calibrationperiod1

global I0 npatches onset_thr

flag1=flagX;

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global method1
global dist1
global factor1

global smoothfactor1

global LBe UBe

close all

% <=============================================================================================>
% <========== Load previously inferred models according to AICc (best to worst fits) ===========>
% <=============================================================================================>

%load(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4),'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

% remove repeated rows
[RMSES,index1]=unique(RMSES,'rows','stable');
PS=PS(index1,:);

index1=rank1;

%[npatches,onset_thr,P0]

% RMSES and PS are sparse; the fit functions work with full values.
npatches=full(RMSES(index1,1));

onset_thr=full(RMSES(index1,2));

AICc_best=full(RMSES(index1,3));

P0=full(PS(index1,1:npatches*4+2));

numparams=get_nparams(method1,dist1,npatches,flag1,1,onset_fixed);

%

% fitting

Phats=[];

goodnessfit=[];

modelfits=[];

SSstats1=[];

SSstats2=[];

%data1=load(strcat('./output/',datafilename1));

data1=data1(epidemic_period,:);

data1(:,2)=data1(:,2);

I0=data1(1,2); % initial condition

if I0==0
    
    data1=data1(2:end,:);
    
end

data=data1(:,2);

[max1,index1]=max(data);

timevect=(data1(:,1));

% <=============================================================================================>
% <===== The ranked fit is bootstrapped as it stands (revision 2026-11) ========================>
% <=============================================================================================>
% Before, the ranked vector P0 was optimized again here from a single start.
% When that run moved, the bootstrap was centred on a fit that no longer had
% the AICc it was ranked with. The screen now returns converged fits
% (fitSubepidemicModel.m), so P0 is used unchanged and its score is checked.

d=P0(end); % legacy output: the power d of the ranked fit

fitSettings=subepidemicFitSettings();
odeOptions=subepidemicOdeOptions();

% Same bounds and coordinate scales as the screen used for this candidate.
[LB,UB,isLog]=subepidemicSearchSpace(npatches,data,flag1,method1,fitSettings);

I0=data(1); % initial condition

hold on

% Same series as the screen fitted (raw counts for count likelihoods).
ydata=subepidemicFitSeries(data,smoothfactor1,method1,fitSettings);

% Explicit fit context: the objective reads nothing from globals.
context=struct('timevect',timevect,'ydata',ydata,'I0',I0, ...
    'npatches',npatches,'onset_fixed',onset_fixed,'onset_thr',onset_thr, ...
    'flag1',flag1,'method1',method1,'odeOptions',odeOptions);

[fvalRanked,rankedDetail]=subepidemicObjective(P0,context);
if ~rankedDetail.success
    error('SubEpiPredict:RankedFitInvalid', ...
        'The ranked model %d could not be reproduced: %s',rank1,rankedDetail.message);
end
if rankedDetail.nActive~=npatches
    error('SubEpiPredict:BestModelActivationMismatch', ...
        'Ranked model %d was fitted with %d sub-epidemics but reconstructs with %d active.', ...
        rank1,npatches,rankedDetail.nActive);
end
AICcCheck=getAICc(method1,dist1,npatches,flag1,1,fvalRanked,length(ydata),onset_fixed);
if abs(AICcCheck-AICc_best)>1e-8*max(1,abs(AICc_best))
    subepidemicWarning('SubEpiPredict:RankedScoreChanged', ...
        ['Ranked model %d: AICc %.10g at the saved vector differs from the ranked AICc %.10g. ', ...
         'Data, smoothing or ODE settings changed between the screen and the bootstrap.'], ...
        rank1,AICcCheck,AICc_best);
end

% A ranked fit with r or K at a search limit is set by the limit, not by the data.
rankedAtBound=subepidemicBoundSummary(P0,LB,UB,isLog,npatches);
for k=1:numel(rankedAtBound)
    if rankedAtBound(k).isSearchLimit
        subepidemicWarning('SubEpiPredict:RankedModelAtBound', ...
            ['Ranked model %d has %s at its %s limit (%g). The estimate is set by the limit. ', ...
             'Set %s in options_fit_settings.m to move it.'], ...
            rank1,rankedAtBound(k).name,rankedAtBound(k).side,rankedAtBound(k).limit,rankedAtBound(k).setting);
    end
end

P=P0;
fval=fvalRanked;

rs_hat=P(1,1:npatches);
ps_hat=P(1,npatches+1:2*npatches);
as_hat=P(1,2*npatches+1:3*npatches);
Ks_hat=P(1,3*npatches+1:4*npatches);

alpha_hat=P(end-1);
d_hat=P(end);

% Shared initialization, activation reset, and incidence conversion.
[~,x,totinc,patchIncidence]=simulateSubepidemic(timevect,I0,npatches,onset_fixed, ...
    onset_thr,flag1,rs_hat,ps_hat,as_hat,Ks_hat,odeOptions);


figure(10)
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
nexttile(1)
for j=1:npatches
    
    incidence1=patchIncidence(:,j);
    
    plot(timevect,incidence1)
    hold on
    
end



bestfit=totinc;


if method1==0 & dist1==0 % normal distribution of the error structure

    var1=sum((bestfit-data).^2)./(length(bestfit)-numparams); % last revised: 01 June 2022
    factor1=sqrt(var1);

elseif method1==6 & dist1==6

    var1=sum((bestfit-data).^2)./(length(bestfit)-numparams); % last revised: 04 Oct 2024
    factor1=sqrt(var1)/sqrt(2);  % scale parameter: b= sigma/sqrt(2);

elseif method1==3 | method1==4 | method1==5

    % Revision 2026-11: the synthetic datasets of this ranked model are drawn
    % with ITS OWN fitted dispersion. Before, factor1 kept the value the
    % screen left in the global, i.e. the alpha of the top-ranked model, for
    % every rank.
    factor1=alpha_hat;

end




        
plot(timevect,totinc,'r')
hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('incidence fit')

title(strcat('Num. Subepidemics=',num2str(npatches),'; AICc=',num2str(AICc_best,6)))

legend(strcat('Sub-epidemics=',num2str(npatches),'; C_{thr}=',num2str(onset_thr)))
set(gca,'FontSize',24)
set(gcf,'color','white')


% <=============================================================================================>
% <======================== Vector with parameter estimates ====================================>
% <=============================================================================================>

Ptrue=[rs_hat ps_hat as_hat Ks_hat alpha_hat d_hat];
fprintf('\nRanked model %d, estimates [r p a K alpha d]: %s\n',rank1,mat2str(full(Ptrue),5));
fprintf('Bootstrap: %d refits of ranked model %d ...\n',M,rank1);

% <=============================================================================================>
% <========== Generate parameter uncertainty via parametric bootstrapping ======================>
% <=============================================================================================>

yi=cumsum(totinc);

z=Ptrue; %r p a Ks alpha(MLE-Neg binomial) d

Phatss=zeros(M,4*npatches+2);
SSEs=zeros(M,1);
curves=zeros(length(yi),M);


% <================================================================================================>
% <=========================== Bootstrap refits ===================================================>
% <================================================================================================>
% Each synthetic dataset is refitted by fitSubepidemicModel.m from the ranked
% vector, one jittered copy of it and fitSettings.bootstrapStarts
% Latin-hypercube points (default 1): three local runs, as before, but in
% the scaled search space. The retention policy is unchanged: a refit with
% no converged local minimum is counted as unresolved and its best valid end
% point is kept. bootstrapDiagnostics records what happened to each refit.

bootstrapDiagnostics=struct('success',false(M,1),'usedFallback',false(M,1), ...
    'exitflag',NaN(M,1),'nActive',NaN(M,1),'startAgreement',NaN(M,1), ...
    'gradientInfNorm',NaN(M,1),'polishRun',false(M,1),'polishImproved',false(M,1), ...
    'nAtBound',NaN(M,1),'fval',NaN(M,1));

bootContext=context;

% Bootstrap counts for user-facing diagnostics only.
acceptedBootstrap = 0;

for real=1:M
     
    yirData=AddPoissonError(yi,1,dist1,factor1,d_hat);
    
    curves(:,real)=yirData;
    
    I0=data(1); % initial condition
    
    ydata=yirData;

    bootContext.ydata=yirData;

    [Pboot,fvalBoot,bootInfo]=fitSubepidemicModel(bootContext,P,LB,UB,isLog,[], ...
        fitSettings.bootstrapStarts,fitSettings);

    if bootInfo.success
        nActiveBoot=bootInfo.nActive;
    elseif ~isempty(bootInfo.fallback)
        % unresolved refit: keep its best valid end point (legacy policy)
        Pboot=bootInfo.fallback.P;
        fvalBoot=bootInfo.fallback.fval;
        nActiveBoot=bootInfo.fallback.nActive;
        bootstrapDiagnostics.usedFallback(real)=true;
    else
        error('SubEpiPredict:BootstrapFitFailed', ...
            'Bootstrap dataset %d of ranked model %d has no valid parameter vector: %s', ...
            real,rank1,bootInfo.message);
    end

    acceptedBootstrap = acceptedBootstrap + double(bootInfo.success);

    bootstrapDiagnostics.success(real)=bootInfo.success;
    bootstrapDiagnostics.exitflag(real)=bootInfo.exitflag;
    bootstrapDiagnostics.nActive(real)=nActiveBoot;
    bootstrapDiagnostics.startAgreement(real)=bootInfo.startAgreement;
    bootstrapDiagnostics.gradientInfNorm(real)=bootInfo.gradientInfNorm;
    bootstrapDiagnostics.polishRun(real)=bootInfo.polish.run;
    bootstrapDiagnostics.polishImproved(real)=bootInfo.polish.improved;
    bootstrapDiagnostics.nAtBound(real)=numel(bootInfo.atLowerBound)+numel(bootInfo.atUpperBound);
    bootstrapDiagnostics.fval(real)=fvalBoot;
    
    % Revision m: the same fixed order of interchangeable sub-epidemics as
    % in the screen, so a refit cannot swap their labels.
    Phatss(real,:)=subepidemicOrderLabels(Pboot,npatches,onset_fixed,onset_thr,data(1));
    
end

% Print bootstrap counts only; do not change fitting or retained results.
fprintf('\nBootstrap summary for ranked model %d:\n', rank1);
fprintf('Requested bootstrap datasets: %d\n', M);
fprintf('Accepted bootstrap datasets: %d\n', acceptedBootstrap);
fprintf('Unresolved bootstrap datasets: %d\n', M - acceptedBootstrap);
fprintf('Refits with fewer than %d active sub-epidemics: %d\n', npatches, sum(bootstrapDiagnostics.nActive<npatches));
% Revision m: one rule for "at a bound", that of the PII table
% (subepidemicBoundHitFractions.m: within 0.1% of the bound).
[boundReport,boundLines,nRefitsAtBound]=subepidemicBoundSummary(Phatss,LB,UB,isLog,npatches);
fprintf('Refits with a parameter at a bound: %d\n', nRefitsAtBound);
for k=1:numel(boundLines)
    fprintf('   %s refits\n',boundLines{k});
end
fprintf('Refits in which the polish ran: %d\n\n', sum(bootstrapDiagnostics.polishRun));

bootstrapDiagnostics.rankedAtBound=rankedAtBound;
bootstrapDiagnostics.boundReport=boundReport;

% <================================================================================================>
% <=========================== Parameter identifiability (PII) ====================================>
% <================================================================================================>
% The measures of the QuantDiffForecast toolbox for every parameter of this
% ranked model (subepidemicIdentifiability.m): PII = log10(UB/LB) of the 95%
% bootstrap interval, the share of the search range the interval covers, the
% largest PII the bounds allow and the ratio to it, the share of the refits
% at each bound, and the Monte Carlo error of the interval and of the PII.
% The search bounds are saved with them, so plotFit_subepidemicFramework.m
% reports against the bounds this fit used.
bootstrapDiagnostics.searchLB=LB;
bootstrapDiagnostics.searchUB=UB;
bootstrapDiagnostics.searchIsLog=isLog;
bootstrapDiagnostics.identifiability=subepidemicIdentifiability(Phatss,LB,UB,npatches,P0,fitSettings.PIIscale);
subepidemicPrintPII(bootstrapDiagnostics.identifiability, ...
    sprintf('Parameter identifiability (PII) of ranked model %d',rank1));

% Revision m. A parameter whose 95% interval covers at least
% notIdentifiedRatio (default 0.5) of what the search bounds allow is not
% identified by the data: its interval is the search range, wherever the
% bounds are put. It gets its own warning, and the advice to move a bound
% (below) is kept for parameters that are identified apart from the bound.
piiReport=bootstrapDiagnostics.identifiability;
notIdentified=find(~piiReport.fixed(:) & piiReport.PII_ratio(:)>=fitSettings.notIdentifiedRatio)';
for k=notIdentified
    if piiReport.PII_scale(k)==1
        coverText=sprintf('%.0f%% of the search range on the log scale (%.1f of %.1f decades)', ...
            100*piiReport.PII_ratio(k),piiReport.PII(k),piiReport.PII_box(k));
    else
        coverText=sprintf('%.0f%% of the search range',100*piiReport.PII_ratio(k));
    end
    subepidemicWarning('SubEpiPredict:NotIdentified', ...
        ['Ranked model %d: %s is not identified by these data. Its 95%% interval, %.4g to %.4g, covers %s, ', ...
         'so it reflects the search bounds, not the data; moving a bound moves the interval. ', ...
         'Do not interpret its estimate (%.4g).'], ...
        rank1,piiReport.name{k},piiReport.CI_LB(k),piiReport.CI_UB(k),coverText,piiReport.estimate(k));
end
bootstrapDiagnostics.notIdentified=piiReport.name(notIdentified);
notIdentifiedColumns=piiReport.column(notIdentified);

% More than 2.5% of the refits at a limit of r or K: the 95% percentile
% interval of that parameter ends at the limit instead of where the data put it.
for k=1:numel(boundReport)
    if boundReport(k).isSearchLimit && boundReport(k).fraction>0.025 && ~any(notIdentifiedColumns==boundReport(k).index)
        subepidemicWarning('SubEpiPredict:BootstrapAtBound', ...
            ['Ranked model %d: %d of %d bootstrap refits have %s at its %s limit (%g), so its 95%% interval ', ...
             'ends at the limit. Set %s in options_fit_settings.m to move it.'], ...
            rank1,boundReport(k).count,boundReport(k).total,boundReport(k).name,boundReport(k).side, ...
            boundReport(k).limit,boundReport(k).setting);
    end
end

