% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function performance=plotFit_subepidemicFramework(outbreakx_pass,caddate1_pass)

% Plot model fits and derive performance metrics during the calibration period for the best fitting models

close all

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1

global smoothfactor1

global calibrationperiod1

% <============================================================================>
% <================== Load the parameter values ===============================>
% <============================================================================>

% options.m
[cumulative1_INP,outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options;

% <============================================================================>
% <================================ Dataset ======================================>
% <============================================================================>


if exist('outbreakx_pass','var')==1 & isempty(outbreakx_pass)==0

    outbreakx=outbreakx_pass;

else
    outbreakx=outbreakx_INP;

end

if exist('caddate1_pass','var')==1 & isempty(caddate1_pass)==0

    caddate1=caddate1_pass;
else
    caddate1=caddate1_INP;

end

cadregion=cadregion_INP; % string indicating the region of the time series (USA, Chile, Mexico, Nepal, etc)

caddisease=caddisease_INP;

datatype=datatype_INP;

datevecfirst1=datevecfirst1_INP;

datevecend1=datevecend1_INP;

DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data, 365=yearly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end

cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


% <============================================================================>
% <============================Adjustments to data =================================>
% <============================================================================>

%smoothfactor1=7; % <smoothfactor1>-day rolling average of the case series

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0; %Type of estimation method: 0 = nonlinear least squares

%dist1=0; % Normnal distribution to model error structure

% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

npatchess2=npatches_fixed;  % maximum number of subepidemics considered in epidemic trajectory fit

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flagss2=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory

% <==============================================================================>
% <======== Number of best fitting models used to generate ensemble model ========================>
% <==============================================================================>

topmodels1=1:topmodelsx_INP;

if npatchess2==1
    topmodels1=1;
end

% Revision i: warnings of this run are collected and listed at its end.
subepidemicWarningLog('start');

% Revision h: use the ranked models the saved fit holds, if options.m asks for
% more (subepidemicAvailableTopModels.m).
topmodels1=subepidemicAvailableTopModels(topmodels1,npatches_fixed,onset_fixed,smoothfactor1, ...
    strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1), ...
    flag1_INP,method1,dist1,calibrationperiod1);

% Random seed, if the setting rngSeed is given (subepidemicSeed.m).
rngInfo=subepidemicSeed('plotFit_subepidemicFramework');


% <=======================================================================================>
% <========== Initialize variables to store results across top-ranked models ===========================>
% <=======================================================================================>

RMSES=[];
MAES=[];
MSES=[];
PIS=[];
MISS=[];
WISS=[];

cc1=1;

quantilescs=[];

param_rs=[];
param_ps=[];
param_as=[];
param_K0s=[];
param_qs=[];
param_alphas=[];
param_ds=[];

AICc_rank1=[];
relativelik_rank1=[];

piiRanked={}; % PII metrics of each ranked model (subepidemicIdentifiability.m)

figure(300+topmodels1(end)+1)

tiledlayout(length(topmodels1),3,'Padding','compact','TileSpacing','compact');

for rank1=topmodels1

    cc2=1;

    npatches_fixed=npatchess2;

    flag1=flagss2(cc2,:);


    % <========================================================================================>
    % <================================ Load model results ====================================>
    % <========================================================================================>

    % Do not retain a previous rank's fit if the current MAT file omits it.
    clear bestfit bestfit_original P0 bootstrapDiagnostics
    load (strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat'))

    % Preserve this rank's saved original-data fit for residual diagnostics.
    % Bootstrap trajectories below are separate from this reference curve.
    if ~exist('bestfit','var') || ~isnumeric(bestfit) || ~isreal(bestfit) || ...
            ~isvector(bestfit) || isempty(bestfit) || ...
            numel(bestfit) ~= size(data1,1) || any(~isfinite(bestfit(:)))
        error('plotFit_subepidemicFramework:InvalidSavedBestFit', ...
            ['Rank %d must contain a finite real bestfit vector with ' ...
             'one value per calibration data row.'], rank1);
    end
    bestfit_original = bestfit(:);

    % The ranked fit to the observed data. The titles and console lines below
    % report it as the estimate; the bootstrap refits give the interval.
    if ~exist('P0','var') || ~isnumeric(P0) || ~isreal(P0) || ...
            numel(P0) ~= 4*npatches+2 || any(~isfinite(P0(:)))
        error('plotFit_subepidemicFramework:InvalidSavedEstimate', ...
            ['Rank %d must contain the fitted parameter vector P0 with ' ...
             '%d finite values. Rerun the fit.'], rank1, 4*npatches+2);
    end
    fittedEstimate = full(double(P0(:)'));

     datevecfirst1=datevecfirst1_INP;


    AICc_rank1=[AICc_rank1;[rank1 AICc_best]];
    relativelik_rank1=[relativelik_rank1;[rank1 relativelik_i(rank1)]];

    timevect=data1(:,1);

    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>


    rs=Phatss(:,1:npatches);
    ps=Phatss(:,npatches+1:2*npatches);
    as=Phatss(:,2*npatches+1:3*npatches);
    Ks=Phatss(:,3*npatches+1:4*npatches);


    % Parameter values
    param_r=[];
    param_p=[];
    param_a=[];
    param_K=[];


    param_alpha=[mean(Phatss(:,end-1)) quantile(Phatss(:,end-1),0.025) quantile(Phatss(:,end-1),0.975)];

    param_d=[mean(Phatss(:,end)) quantile(Phatss(:,end),0.025) quantile(Phatss(:,end),0.975)];

    % Estimate of the ranked fit with the bootstrap interval (param_alpha and
    % param_d keep the bootstrap mean for the csv files below).
    cad6=subepidemicParameterLabel('alpha',fittedEstimate(end-1),param_alpha(end,2),param_alpha(end,3))

    cad7=subepidemicParameterLabel('d',fittedEstimate(end),param_d(end,2),param_d(end,3))

    i=1;

    figure(300+rank1)

    tiledlayout(full(double(npatches)),4,'Padding','compact','TileSpacing','compact');

    for j=1:npatches

        param_r=[param_r;[mean(rs(:,j)) quantile(rs(:,j),0.025) quantile(rs(:,j),0.975)]];

        param_p=[param_p;[mean(ps(:,j)) quantile(ps(:,j),0.025) quantile(ps(:,j),0.975)]];

        param_a=[param_a;[mean(as(:,j)) quantile(as(:,j),0.025) quantile(as(:,j),0.975)]];

        param_K=[param_K;[mean(Ks(:,j)) quantile(Ks(:,j),0.025) quantile(Ks(:,j),0.975)]];

        % Titles: the estimate of the ranked fit (not the bootstrap mean, which
        % differs from it for a skewed parameter) and the bootstrap interval,
        % three significant digits. param_r ... param_K keep the bootstrap mean
        % in their first column; the csv files below report mean, estimate and
        % median.
        cad1=subepidemicParameterLabel(sprintf('r_%d',j),fittedEstimate(j),param_r(j,2),param_r(j,3));
        cad2=subepidemicParameterLabel(sprintf('p_%d',j),fittedEstimate(npatches+j),param_p(j,2),param_p(j,3))
        cad3=subepidemicParameterLabel(sprintf('a_%d',j),fittedEstimate(2*npatches+j),param_a(j,2),param_a(j,3))
        cad4=subepidemicParameterLabel(sprintf('K_%d',j),fittedEstimate(3*npatches+j),param_K(j,2),param_K(j,3))


        nexttile(i)

        hist(rs(:,j))
        hold on

        xlabel('r')
        title(cad1)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_r(j,2) 10;param_r(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--');
        set(line1,'LineWidth',2)

        % fitted estimate of r_j: black dashed vertical line
        subepidemicEstimateLine(fittedEstimate(j));

        legend(strcat(num2ordinal(rank1),' Ranked Model'))

        nexttile(i+1)
        hist(ps(:,j))
        hold on

        xlabel('p')
        title(cad2)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_p(j,2) 10;param_p(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--');
        set(line1,'LineWidth',2)

        % fitted estimate of p_j: black dashed vertical line
        subepidemicEstimateLine(fittedEstimate(npatches+j));


        nexttile(i+2)
        hist(as(:,j))
        hold on
        xlabel('a')
        title(cad3)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_a(j,2) 10;param_a(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--');
        set(line1,'LineWidth',2)

        % fitted estimate of a_j: black dashed vertical line
        subepidemicEstimateLine(fittedEstimate(2*npatches+j));


        nexttile(i+3)
        hist(Ks(:,j))
        hold on
        xlabel('K')
        title(cad4)

        line2=[param_K(j,2) 10;param_K(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--');
        set(line1,'LineWidth',2)

        % fitted estimate of K_j: black dashed vertical line
        subepidemicEstimateLine(fittedEstimate(3*npatches+j));


        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        i=i+4;


    end

    % <==============================================================================================>
    % <=============== Save csv files with parameter estimates for each sub-epidemic ================>
    % <==============================================================================================>

    forecastingperiod=0;

    % Columns 1 to 5 are as before: "mean" is the mean of the bootstrap refits.
    % SCI is the log10 width of the interval, computed by sci_interval.m, which
    % returns 6 instead of Inf when the lower limit is 0.
    % Two columns are added at the end: "estimate", the value in the ranked fit
    % to the observed data (the number in the histogram title), and "median",
    % the median of the refits. A few refits with a very large value can move
    % the mean far from both.

    T = array2table([(1:npatches)' param_r sci_interval(param_r(:,2),param_r(:,3)) fittedEstimate(1:npatches)' median(Phatss(:,1:npatches),1)']);
    T.Properties.VariableNames(1)={'sub-epidemic'};
    T.Properties.VariableNames(2:7) = {'r mean','r 95% CI LB','r 95% CI UB','r SCI','r estimate','r median'};
    writetable(T,strcat('./output/param-r-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    T = array2table([(1:npatches)' param_p sci_interval(param_p(:,2),param_p(:,3)) fittedEstimate(npatches+1:2*npatches)' median(Phatss(:,npatches+1:2*npatches),1)']);
    T.Properties.VariableNames(1)={'sub-epidemic'};
    T.Properties.VariableNames(2:7) = {'p mean','p 95% CI LB','p 95% CI UB','p SCI','p estimate','p median'};
    writetable(T,strcat('./output/param-p-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    T = array2table([(1:npatches)' param_a sci_interval(param_a(:,2),param_a(:,3)) fittedEstimate(2*npatches+1:3*npatches)' median(Phatss(:,2*npatches+1:3*npatches),1)']);
    T.Properties.VariableNames(1)={'sub-epidemic'};
    T.Properties.VariableNames(2:7) = {'a mean','a 95% CI LB','a 95% CI UB','a SCI','a estimate','a median'};
    writetable(T,strcat('./output/param-a-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    T = array2table([(1:npatches)' param_K sci_interval(param_K(:,2),param_K(:,3)) fittedEstimate(3*npatches+1:4*npatches)' median(Phatss(:,3*npatches+1:4*npatches),1)']);
    T.Properties.VariableNames(1)={'sub-epidemic'};
    T.Properties.VariableNames(2:7) = {'K mean','K 95% CI LB','K 95% CI UB','K SCI','K estimate','K median'};
    writetable(T,strcat('./output/param-K-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

    % <==============================================================================================>
    % <=============== Save csv files with N.B error structure parameters ===========================>
    % <==============================================================================================>

    if method1==3 || method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;

        T = array2table([param_alpha sci_interval(param_alpha(:,2),param_alpha(:,3)) fittedEstimate(end-1) median(Phatss(:,end-1),1)]);
        T.Properties.VariableNames(1:6) = {'alpha mean','alpha 95% CI LB','alpha 95% CI UB','alpha SCI','alpha estimate','alpha median'};
        writetable(T,strcat('./output/param-NB-alpha-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))
    elseif method1==5
        T = array2table([param_alpha sci_interval(param_alpha(:,2),param_alpha(:,3)) fittedEstimate(end-1) median(Phatss(:,end-1),1)]);
        T.Properties.VariableNames(1:6) = {'alpha mean','alpha 95% CI LB','alpha 95% CI UB','alpha SCI','alpha estimate','alpha median'};
        writetable(T,strcat('./output/param-NB-alpha-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))

        T = array2table([param_d sci_interval(param_d(:,2),param_d(:,3)) fittedEstimate(end) median(Phatss(:,end),1)]);
        T.Properties.VariableNames(1:6) = {'d mean','d 95% CI LB','d 95% CI UB','d SCI','d estimate','d median'};
        writetable(T,strcat('./output/param-NB-d-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))
    end


    % <==============================================================================================>
    % <=============== Parameter identifiability (PII) of this ranked model =========================>
    % <==============================================================================================>
    % The measures the QuantDiffForecast toolbox reports (subepidemicIdentifiability.m),
    % one row per parameter:
    %   PII            log10(UB/LB) of the 95% bootstrap interval; NaN when the lower
    %                  limit is <= 0 (the SCI column above gives 6 there)
    %   PII_rangeFrac  interval width as a share of the search range, in [0,1]
    %   PII_scale      headline measure: 1 = PII (lower bound > 0), 0 = PII_rangeFrac
    %   PII_box        log10(UB/LB) of the search bounds, the largest PII they allow
    %   PII_ratio      ~1 => the interval is the search range (censored by the bounds)
    %   fracLB, fracUB share of the bootstrap refits at the lower, upper bound
    %   MCSE_*         Monte Carlo SE of the interval limits and of the PII
    % The bounds are those saved by the fit. A fit saved before this report
    % existed does not hold them; they are then rebuilt from the current settings.
    piiSettings=subepidemicFitSettings();
    if exist('bootstrapDiagnostics','var') && isstruct(bootstrapDiagnostics) && ...
            isfield(bootstrapDiagnostics,'searchLB') && isfield(bootstrapDiagnostics,'searchUB') && ...
            numel(bootstrapDiagnostics.searchLB)==4*npatches+2 && numel(bootstrapDiagnostics.searchUB)==4*npatches+2
        piiLB=bootstrapDiagnostics.searchLB;
        piiUB=bootstrapDiagnostics.searchUB;
    else
        [piiLB,piiUB]=subepidemicSearchSpace(npatches,data1(:,2),flag1(1),method1,piiSettings);
        subepidemicWarning('SubEpiPredict:PIIBoundsRebuilt', ...
            ['Ranked model %d was saved without its search bounds. The PII report uses the bounds of the ', ...
             'current settings; if the bound settings changed since the fit, rerun the fit.'],rank1);
    end

    piiRanked{end+1}=subepidemicIdentifiability(Phatss,piiLB,piiUB,npatches,fittedEstimate,piiSettings.PIIscale); %#ok<AGROW>

    subepidemicPrintPII(piiRanked{end},sprintf('Parameter identifiability (PII) of ranked model %d',rank1));

    subepidemicWritePII(strcat('./output/PII-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'),piiRanked{end});


    % <========================================================================================>
    % <================================ Plot model fit ========================================>
    % <========================================================================================>

    LB1=quantile(curves',0.025);
    LB1=(LB1>=0).*LB1;

    UB1=quantile(curves',0.975);
    UB1=(UB1>=0).*UB1;

    figure(300+topmodels1(end)+1)

    nexttile(cc1)

    plot(timevect,curves,'c-')
    hold on

    line1=plot(timevect,median(curves,2),'r-');

    set(line1,'LineWidth',2)

    line1=plot(timevect,LB1,'r--');
    set(line1,'LineWidth',2)

    line1=plot(timevect,UB1,'r--');
    set(line1,'LineWidth',2)

    line1=plot(timevect,data1(:,2),'ko');
    set(line1,'LineWidth',2)


    hold on

    xlabel('Time')
    ylabel(strcat(caddisease,{' '},datatype))

    axis([timevect(1) timevect(end)+1 0 max(data(:,2))*1.3])

    set(gca,'FontSize', 16);
    set(gcf,'color','white')




    % <=============================================================================================>
    % <================================ Plot subepidemic profile ===========================================>
    % <=============================================================================================>

    nexttile(cc1+1)

    color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];

    fittedcurves=zeros(length(timevect),M);

    % generate forecast curves from each bootstrap realization
    for realization=1:M

        rs_hat=Phatss(realization,1:npatches);
        ps_hat=Phatss(realization,npatches+1:2*npatches);
        as_hat=Phatss(realization,2*npatches+1:npatches*3);
        Ks_hat=Phatss(realization,3*npatches+1:4*npatches);

        % Shared initialization, activation reset, and incidence conversion.
        [~,x,totinc,patchIncidence]=simulateSubepidemic(timevect,data1(1,2),npatches,onset_fixed, ...
            onset_thr,flag1,rs_hat,ps_hat,as_hat,Ks_hat);


        for j=1:npatches

            incidence1=patchIncidence(:,j);

            plot(timevect,incidence1,color1(j,:))

            hold on

        end



        % Keep bootstrap curves separate from the saved original-data fit.

        fittedcurves(:,realization)=totinc;

        gray1=gray(10);

        plot(timevect,totinc,'color',gray1(7,:))

        %hold on
        %plot(timevect,data1(:,2),'ko')

    end

    xlabel('Time (days)');
    ylabel(strcat(caddisease,{' '},datatype))

    line1=plot(data(:,1),data(:,2),'ko');
    set(line1,'LineWidth',2)

    axis([timevect(1) timevect(end)+1 0 max(data(:,2))*1.3])

    set(gca,'FontSize',16)
    set(gcf,'color','white')

    %title('Sub-epidemic profile')

    %title(strcat('Num. Subepidemics=',num2str(npatches),'; AICc=',num2str(AICc_best,6)))

    legend(strcat('Sub-epidemics=',num2str(npatches),'; C_{thr}=',num2str(onset_thr)))

    title(strcat(num2ordinal(rank1),' Ranked Model'))

    % <========================================================================================>
    % <================================ Store model fit quantiles ======================================>
    % <========================================================================================>

    [quantilesc,quantilesf]=computeQuantiles(data1,curves,0);
    quantilescs=[quantilescs; quantilesc];


    % compute doubling times

    forecastingperiod=0;
    
    % Revision 2026-11-16: meandoublingtime starts empty. It used to start as M zeros
    % and the M values were appended after them, which halved the mean and put the
    % lower limit at 0. Curves with no doubling give NaN and are left out.
    meandoublingtime=[];

    doublingtimess=zeros(30,M)+NaN;

    maxd=1;

    for j=1:M

        [tds,C0data,curve,doublingtimes]=getDoublingTimeCurve(max(fittedcurves(:,j),0),DT,0);

        doublingtimess(1:length(doublingtimes),j)=doublingtimes;

        if maxd<length(doublingtimes)
            maxd=length(doublingtimes);
        end

        meandoublingtime=[meandoublingtime;mean(doublingtimes)];
    end

    doublingtimess=doublingtimess(1:maxd,1:M);

    seq_doublingtimes=[];

    for j=1:maxd

        index1=find(~isnan(doublingtimess(j,:)));

        seq_doublingtimes=[seq_doublingtimes;[j mean(doublingtimess(j,index1)) quantile(doublingtimess(j,index1),0.025) quantile(doublingtimess(j,index1),0.975) length(index1)./M]];

    end

    seq_doublingtimes % [ith doubling, mean, 95%CI LB, 95%CI UB, prob. i_th doubling]

    % Mean doubling times
    meandoublingtime=meandoublingtime(~isnan(meandoublingtime));
    if isempty(meandoublingtime) % no curve doubled
        dmean=NaN;
        dLB=NaN;
        dUB=NaN;
    else
        dmean=mean(meandoublingtime);
        dLB=quantile(meandoublingtime,0.025);
        dUB=quantile(meandoublingtime,0.975);
    end

    param_doubling=[dmean dLB dUB]

    % <=============================================================================================>
    % <============================== Save file with doubling time estimates =======================>
    % <=============================================================================================>

    T = array2table(seq_doublingtimes);
    T.Properties.VariableNames(1:5) = {'i_th doubling','db mean','db 95%CI LB','db 95% CI UB','prob. i_th doubling'};
    writetable(T,strcat('./output/doublingTimes-ranked(', num2str(rank1),')-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))



    % <========================================================================================>
    % <================================ Plot residuals ========================================>
    % <========================================================================================>

    nexttile(cc1+2)

    % Retain the existing fitted-minus-observed sign convention.
    resid1 = bestfit_original - data1(:,2);

    stem(timevect,resid1,'b')
    hold on

    xlabel('Time (days)')
    ylabel('Residuals (fitted - observed)')

    axis([timevect(1) timevect(end)+1 min(resid1)-1 max(resid1)+1])

    set(gca,'FontSize', 16);
    set(gcf,'color','white')


    figure(400)
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');

    hist(Phatss(:,3*npatches+1:4*npatches))
    xlabel('Sub-epidemic size')
    ylabel('Frequency')

    set(gca,'FontSize', 16);
    set(gcf,'color','white')


    mean(Phatss(:,3*npatches+1:4*npatches))

    sum(Phatss(:,3*npatches+1:4*npatches),2);

    mean(sum(Phatss(:,3*npatches+1:4*npatches),2))


    % <========================================================================================>
    % <================================ Get performance metrics ===================================>
    % <============================================================================================>

    [RMSEC MSEC MAEC PIC MISC RMSEF MSEF MAEF PIF MISF]=computeforecastperformance(data1,data1,fittedcurves,curves,0);

    [WISC,WISFS]=computeWIS(data1,data1,curves,0);

    % store calibration performance metrics
    RMSES=[RMSES;RMSEC];
    MAES=[MAES;MAEC];
    MSES=[MSES;MSEC];
    PIS=[PIS;PIC];
    MISS=[MISS;MISC];
    WISS=[WISS;WISC];

    cc2=cc2+1;

    cc1=cc1+3;

end

% <============================================================================>
% <=================plot calibration performance metrics for the top-ranked models ==============>
% <============================================================================>

figure(400)

tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile(1)
line1=plot(MAES,'k-o');
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('MAE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

nexttile(2)
line1=plot(MSES,'k-o');
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('MSE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

nexttile(3)
line1=plot(PIS,'k-o');
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('Coverage of the 95% PI')

set(gca,'FontSize', 16);
set(gcf,'color','white')

nexttile(4)
line1=plot(WISS,'k-o');
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('WIS')

set(gca,'FontSize', 16);
set(gcf,'color','white')

% <=============================================================================================>
% <============================== Save file with top-ranked performance metrics (calibration) ====================>
% <=============================================================================================>


performance=[topmodels1' MAES MSES PIS WISS AICc_rank1(:,2) relativelik_rank1(:,2)];

T = array2table(performance);
T.Properties.VariableNames(1:7) = {'i_th-ranked model','MAE','MSE','Coverage 95%PI','WIS','AICc','RelativeLikelihood'};
writetable(T,strcat('./output/performance-calibration-topRanked-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


% <=============================================================================================>
% <============================== Save file with the PII metrics of all top-ranked models ======>
% <=============================================================================================>
% The rows of the PII-ranked(k) files together, with the rank and the number
% of sub-epidemics of each model in the first two columns.

subepidemicWritePII(strcat('./output/PII-topRanked-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'),piiRanked,topmodels1);


%%
% <==================================================================================>
% <========== compute ensemble model from top fitting models and evaluate performance metrics ============>
% <==================================================================================>

if 0

    if length(topmodels1)>1

        figure(500)

        tiledlayout(1,1,'Padding','compact','TileSpacing','compact');

        forecastingperiod=0;
        getperformance=1;

        weights1=subepidemicEnsembleWeights(1,AICc_bests,[],[]); % Akaike weights; (1/AICc)/sum(1/AICc) is not a weight

        curvesforecasts1ens=[];
        curvesforecasts2ens=[];


        for rank1=topmodels1

            % <========================================================================================>
            % <================================ Load model results ====================================>
            % <========================================================================================>

            load (strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat'))

            M1=M;

            curvesforecasts1=zeros(length(timevect),M1);

            curvesforecasts2=[];

            for j=1:M1

                P=Phatss(j,:);

                rs_hat=P(1,1:npatches);
                ps_hat=P(1,npatches+1:2*npatches);
                as_hat=P(1,2*npatches+1:3*npatches);
                Ks_hat=P(1,3*npatches+1:4*npatches);

                alpha_hat=P(1,end-1);
                d_hat=P(1,end);

                % Shared initialization, activation reset, and incidence conversion.
                [~,x,totinc,patchIncidence]=simulateSubepidemic(timevect,data1(1,2),npatches,onset_fixed, ...
                    onset_thr,flag1,rs_hat,ps_hat,as_hat,Ks_hat);




                %

                fittedCurve1=totinc;

                curvesforecasts1(:,j)= fittedCurve1;

                curvesforecasts2=[curvesforecasts2 AddPoissonError(cumsum(curvesforecasts1(:,j)),20,dist1,factor1,1)];


            end

            M1=length(curvesforecasts1(1,:));

            index1=datasample(1:M1,round(M1*weights1(rank1)),'Replace',false);

            curvesforecasts1ens=[curvesforecasts1ens curvesforecasts1(:,index1)];


            M2=length(curvesforecasts2(1,:));

            index2=datasample(1:M2,round(M2*weights1(rank1)),'Replace',false);

            curvesforecasts2ens=[curvesforecasts2ens curvesforecasts2(:,index2)];


        end

        curvesforecasts1=curvesforecasts1ens;

        curvesforecasts2=curvesforecasts2ens;


        timevect=data1(:,1);


        % <========================================================================================>
        % <================================ Plot ensemble model fit ====================================>
        % <========================================================================================>

        caddate1=datestr(datenum(caddate1),'mm-dd-yyyy');

        datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

        datevec1=datevec(datenum1+forecastingperiod*DT);

        wave=[datevecfirst1 datevec1(1:3)];

        hold on

        quantile(curvesforecasts2',0.025)

        LB1=quantile(curvesforecasts2',0.025);
        UB1=quantile(curvesforecasts2',0.975);

        size(LB1)
        size(timevect)

        h=area(timevect',[LB1' UB1'-LB1'])
        hold on

        h(1).FaceColor = [1 1 1];
        h(2).FaceColor = [0.8 0.8 0.8];

        line1=plot(timevect,median(curvesforecasts2,2),'r-');

        set(line1,'LineWidth',2)

        line1=plot(timevect,LB1,'r--');
        set(line1,'LineWidth',2)

        line1=plot(timevect,UB1,'r--');
        set(line1,'LineWidth',2)


        % plot time series datalatest
        line1=plot(data1(:,1),data1(:,2),'ko');
        set(line1,'LineWidth',2)


        axis([0 length(timevect)-1 0 max(UB1)*1.2])

        line2=[timevect(end) 0;timevect(end) max(UB1)*1.20];

        box on


        % plot dates in x axis
        'day='
        datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
        datestr(datenum1)

        datenumIni=datenum1;
        datenumEnd=datenum(wave(4:6))

        dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');


        if DT==1

            set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

        elseif DT==7

            set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

        elseif DT==365

            years1=wave(1)+timelags:wave(1)+timelags+length(dates1(:,1))-1;

            set(gca,'XTick',0:1:length(years1)-1);

            set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

        end


        xticklabel_rotate;


        ylabel(strcat(caddisease,{' '},datatype))

        title('Ensemble Model')

        set(gca,'FontSize',24)
        set(gcf,'color','white')

        % <==========================================================================================>
        % <==================== Get calibration performance metrics for ensemble model ==============>
        % <==========================================================================================>

        if getperformance

            [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1, RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,data1,curvesforecasts1,curvesforecasts2,forecastingperiod);

            [WISC,WISFS]=computeWIS(data1,data1,curvesforecasts2,forecastingperiod)

            % <=============================================================================================>
            % <============================== Save file with ensemble performance metrics ==============================>
            % <=============================================================================================>

            performanceEns=[topmodels1(end) MAECS_model1  MSEFS_model1 PIFS_model1 WISC];

            T = array2table(performanceEns);
            T.Properties.VariableNames(1:5) = {'Ensemble(i) model','MAE','MSE','Coverage 95%PI','WIS'};
            writetable(T,strcat('./output/performance-calibration-Ensemble-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-horizon-',num2str(forecastingperiod),'-weight_type-',num2str(weight_type1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


        end

    end

end

% Revision i: list the toolbox's warnings of this run once more (subepidemicWarningLog.m).
subepidemicWarningLog('summary','plotFit_subepidemicFramework');
