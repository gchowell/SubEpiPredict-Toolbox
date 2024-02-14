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
[cumulative1_INP,outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options

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

for rank1=topmodels1

    cc2=1;

    npatches_fixed=npatchess2

    flag1=flagss2(cc2,:);


    % <========================================================================================>
    % <================================ Load model results ====================================>
    % <========================================================================================>

    load (strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat'))

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

    cad6=strcat('alpha=',num2str(param_alpha(end,1),3),'(95%CI:',num2str(param_alpha(end,2),3),',',num2str(param_alpha(end,3),3),')')

    cad7=strcat('d=',num2str(param_d(end,1),3),'(95%CI:',num2str(param_d(end,2),3),',',num2str(param_d(end,3),3),')')

    i=1;

    figure(300+rank1)

    for j=1:npatches

        param_r=[param_r;[mean(rs(:,j)) quantile(rs(:,j),0.025) quantile(rs(:,j),0.975)]];

        param_p=[param_p;[mean(ps(:,j)) quantile(ps(:,j),0.025) quantile(ps(:,j),0.975)]];

        param_a=[param_a;[mean(as(:,j)) quantile(as(:,j),0.025) quantile(as(:,j),0.975)]];

        param_K=[param_K;[mean(Ks(:,j)) quantile(Ks(:,j),0.025) quantile(Ks(:,j),0.975)]];

        cad1=strcat('r_',num2str(j),'=',num2str(param_r(j,1),2),'(95%CI:',num2str(param_r(j,2),2),',',num2str(param_r(j,3),2),')');
        cad2=strcat('p_',num2str(j),'=',num2str(param_p(j,1),2),'(95%CI:',num2str(param_p(j,2),2),',',num2str(param_p(j,3),2),')')
        cad3=strcat('a_',num2str(j),'=',num2str(param_a(j,1),2),'(95%CI:',num2str(param_a(j,2),2),',',num2str(param_a(j,3),2),')')
        cad4=strcat('K_',num2str(j),'=',num2str(param_K(j,1),3),'(95%CI:',num2str(param_K(j,2),3),',',num2str(param_K(j,3),3),')')


        subplot(npatches,4,i)

        hist(rs(:,j))
        hold on

        xlabel('r')
        title(cad1)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_r(j,2) 10;param_r(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        legend(strcat(num2ordinal(rank1),' Ranked Model'))

        subplot(npatches,4,i+1)
        hist(ps(:,j))
        hold on

        xlabel('p')
        title(cad2)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_p(j,2) 10;param_p(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)


        subplot(npatches,4,i+2)
        hist(as(:,j))
        hold on
        xlabel('a')
        title(cad3)

        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        line2=[param_a(j,2) 10;param_a(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)


        subplot(npatches,4,i+3)
        hist(Ks(:,j))
        hold on
        xlabel('K')
        title(cad4)

        line2=[param_K(j,2) 10;param_K(j,3) 10];

        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)


        set(gca,'FontSize', 16);
        set(gcf,'color','white')

        i=i+4;

    end

    % <========================================================================================>
    % <================================ Plot model fit ========================================>
    % <========================================================================================>

    LB1=quantile(curves',0.025);
    LB1=(LB1>=0).*LB1;

    UB1=quantile(curves',0.975);
    UB1=(UB1>=0).*UB1;

    figure(300+topmodels1(end)+1)

    subplot(length(topmodels1),3,cc1)

    plot(timevect,curves,'c-')
    hold on

    line1=plot(timevect,median(curves,2),'r-') 

    set(line1,'LineWidth',2)

    line1=plot(timevect,LB1,'r--')
    set(line1,'LineWidth',2)

    line1=plot(timevect,UB1,'r--')
    set(line1,'LineWidth',2)

    line1=plot(timevect,data1(:,2),'ko')
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

    subplot(length(topmodels1),3,cc1+1)

    color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];

    fittedcurves=zeros(length(timevect),M);

    % generate forecast curves from each bootstrap realization
    for realization=1:M

        rs_hat=Phatss(realization,1:npatches);
        ps_hat=Phatss(realization,npatches+1:2*npatches);
        as_hat=Phatss(realization,2*npatches+1:npatches*3);
        Ks_hat=Phatss(realization,3*npatches+1:4*npatches);

        IC=zeros(npatches,1);

        if onset_thr>0
            IC(1,1)=data1(1,2);
            IC(2:end,1)=1;

            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);

            invasions(1)=1;
            timeinvasions(1)=0;
            Cinvasions(1)=0;
        else
            IC(1:end,1)=data1(1,2)./length(IC(1:end,1));

            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);

            invasions(1:end)=1;
            timeinvasions(1:end)=0;
            Cinvasions(1:end)=0;
        end


        [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);


        for j=1:npatches

            incidence1=[x(1,j);diff(x(:,j))];

            plot(timevect,incidence1,color1(j,:))

            hold on

        end

        y=sum(x,2);

        totinc=[y(1,1);diff(y(:,1))];

        totinc(1)=totinc(1)-(npatches-1);

        bestfit=totinc;

        fittedcurves(:,realization)=totinc;

        gray1=gray(10);

        plot(timevect,totinc,'color',gray1(7,:))

        %hold on
        %plot(timevect,data1(:,2),'ko')

    end

    xlabel('Time (days)');
    ylabel(strcat(caddisease,{' '},datatype))

    line1=plot(data(:,1),data(:,2),'ko')
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
    
    meandoublingtime=zeros(M,1);

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
    dmean=mean(meandoublingtime);
    dLB=quantile(meandoublingtime,0.025);
    dUB=quantile(meandoublingtime,0.975);

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

    subplot(length(topmodels1),3,cc1+2)

    resid1=bestfit-data1(:,2);

    stem(timevect,resid1,'b')
    hold on

    xlabel('Time (days)')
    ylabel('Residuals')

    axis([timevect(1) timevect(end)+1 min(resid1)-1 max(resid1)+1])

    set(gca,'FontSize', 16);
    set(gcf,'color','white')


    figure(400)
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

subplot(2,2,1)
line1=plot(MAES,'k-o')
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('MAE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,2)
line1=plot(MSES,'k-o')
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('MSE')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,3)
line1=plot(PIS,'k-o')
set(line1,'linewidth',2)
xlabel('i_{th}Ranked Model')
ylabel('Coverage of the 95% PI')

set(gca,'FontSize', 16);
set(gcf,'color','white')

subplot(2,2,4)

line1=plot(WISS,'k-o')
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


%%
% <==================================================================================>
% <========== compute ensemble model from top fitting models and evaluate performance metrics ============>
% <==================================================================================>

if 0

    if length(topmodels1)>1

        figure(500)

        forecastingperiod=0;
        getperformance=1;

        weights1=(1./AICc_bests)/(sum(1./AICc_bests));

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

                IC=zeros(npatches,1);

                if onset_thr>0
                    IC(1,1)=data1(1,2);
                    IC(2:end,1)=1;

                    invasions=zeros(npatches,1);
                    timeinvasions=zeros(npatches,1);
                    Cinvasions=zeros(npatches,1);

                    invasions(1)=1;
                    timeinvasions(1)=0;
                    Cinvasions(1)=0;
                else
                    IC(1:end,1)=data1(1,2)./length(IC(1:end,1));

                    invasions=zeros(npatches,1);
                    timeinvasions=zeros(npatches,1);
                    Cinvasions=zeros(npatches,1);

                    invasions(1:end)=1;
                    timeinvasions(1:end)=0;
                    Cinvasions(1:end)=0;
                end


                [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);



                y=sum(x,2);

                totinc=[y(1,1);diff(y(:,1))];

                if onset_fixed==0
                    totinc(1)=totinc(1)-(npatches-1);
                end
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

        line1=plot(timevect,median(curvesforecasts2,2),'r-')

        set(line1,'LineWidth',2)

        line1=plot(timevect,LB1,'r--')
        set(line1,'LineWidth',2)

        line1=plot(timevect,UB1,'r--')
        set(line1,'LineWidth',2)


        % plot time series datalatest
        line1=plot(data1(:,1),data1(:,2),'ko')
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

        end

    end

end

