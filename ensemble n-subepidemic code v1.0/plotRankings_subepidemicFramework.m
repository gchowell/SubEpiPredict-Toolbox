
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function plotRankings_subepidemicFramework(outbreakx_pass,caddate1_pass)

% Plot ranked models and their AICc related values


% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1

global smoothfactor1 calibrationperiod1

% <============================================================================>
% <================== Load the parameter values ===============================>
% <============================================================================>

% options.m
[cumulative1_INP, outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options

% <============================================================================>
% <================================ Dataset ====================================>
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

cadregion=cadregion_INP;

caddisease=caddisease_INP;

datatype=datatype_INP;

DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end


cadfilename1=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

%smoothfactor1=7; % <smoothfactor1>-day smoothing

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0;

%dist1=0;

% <==============================================================================>
% <========================= Mathematical model =================================>
% <==============================================================================>

%npatches_fixed=2;

%type of growth model used to describe a subepidemic

% 0 = GGM
% 1 = GLM
% 2 = GRM
% 3 = LM
% 4 = Richards

flag1=flag1_INP; % Sequence of subepidemic growth models considered in epidemic wave

% <==============================================================================>
% <======== Number of best fitting models used to generate ensemble model =======>
% <==============================================================================>

topmodelsx=1:topmodelsx_INP;

AICc_bests=[];

for run_id=-1
    %for run_id=0:1:58

    cc1=1;

    close all

    %i=(run_id)*30+1;
    %ARIMA_mean1=ARIMAforecasts(i:1:i+forecastingperiod-1,1);
    %ARIMA_lb1=ARIMAforecasts(i:1:i+forecastingperiod-1,11);
    %ARIMA_ub1=ARIMAforecasts(i:1:i+forecastingperiod-1,end-1);


    %outbreakx=52;

    %caddate1='04-20-20';
    %caddate1='05-11-20'; % for paper practical use

    %caddate1='06-29-20';
    %caddate1='05-11-20';
    %caddate1='03-22-21';

    %caddate1='11-15-21';

    cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


    %

    load (strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'))

    
    % remove repeated rows
    [RMSES,index1]=unique(RMSES,'rows','stable');
    PS=PS(index1,:);
    
    
    AICc_bests=[AICc_bests full(RMSES(topmodelsx,3))];
    
    
    if 1  %plot relative likelihoods of the models
        
        figure(99)
        
        subplot(1,3,1)
        line1=plot(RMSES(topmodelsx,3),'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} Ranked Model')
        ylabel('AICc')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        subplot(1,3,2)
        
        line1=plot(relativelik_i(topmodelsx),'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} Ranked Model')
        ylabel('Relative likelihood')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        subplot(1,3,3)
        
        %deltas=AICc_bests-AICc_bests(1);
        
        line1=plot(1./relativelik_i(topmodelsx),'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} Ranked Model')
        ylabel('Evidence ratio')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
    end
    
    
    figure(100)
    
    color1=['r-';'b-';'g-';'m-';'c-';'y-';'r-';'b-';'g-';'m-';'c-';'y-';'r-';'b-';'g-';'m-';'c-';'y-';'r-';'b-';'g-';'m-';'c-';'y-';'r-';'b-';'g-';'m-';'c-';'y-';];
    
    
    factors=factor(length(topmodelsx));
    
    if length(factors)==1
        rows=factors;
        cols=1;
        
    elseif length(factors)==3
        rows=factors(1)*factors(2);
        cols=factors(3);
    else
        rows=factors(1);
        cols=factors(2);
    end
    
    
    for index1=topmodelsx
        
   
        subplot(rows,cols,index1)
        
        npatches=RMSES(index1,1);
        
        onset_thr=RMSES(index1,2);
        
        AICc_best=RMSES(index1,3);
        
        
        
        P=PS(index1,1:npatches*4+2);
        
        rs_hat=P(1,1:npatches)
        ps_hat=P(1,npatches+1:2*npatches)
        as_hat=P(1,2*npatches+1:3*npatches)
        Ks_hat=P(1,3*npatches+1:4*npatches)
        
        alpha_hat=P(1,end-1);
        d_hat=P(1,end);
        
        
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
        
        
        
        IC=zeros(npatches,1);
        
        if onset_fixed==0
            IC(1,1)=I0;
            IC(2:end,1)=1;
            
            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);
            
            invasions(1)=1;
            timeinvasions(1)=0;
            Cinvasions(1)=0;
            
        else
            
            IC(1:end,1)=I0./length(IC(1:end,1));
            
            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);
            
            invasions(1:end)=1;
            timeinvasions(1:end)=0;
            Cinvasions(1:end)=0;
        end
        
        
        [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);
        %x=ode5(@modifiedLogisticGrowthPatch,timevect,IC,r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1);
        

        
        %     if sum(invasions)<npatches
        %
        %         invasions
        %
        %
        %         npatches=sum(invasions);
        %
        %         P=[rs_hat(1:npatches) ps_hat(1:npatches) as_hat(1:npatches) Ks_hat(1:npatches) alpha_hat d_hat];
        %
        %         PS(index1,:)=0;
        %
        %         PS(index1,1:length(P))=P;
        %
        %         RMSES(index1,1)=npatches;
        %
        %     end
        %
        
        for j=1:npatches
            
            incidence1=[x(1,j);diff(x(:,j))];
            
            line1=plot(timevect,incidence1,color1(j,:))
            
            set(line1,'Linewidth',2)
            hold on
            
        end
        
        y=sum(x,2);
        
        totinc=[y(1,1);diff(y(:,1))];
        
        if onset_thr>0
            totinc(1)=totinc(1)-(npatches-1);
        end
        
        bestfit=totinc;
        
        
        hold on
        line1=plot(timevect,data,'bo')
        set(line1,'Linewidth',2,'markerSize',8)
        
        
        line1=plot(timevect,totinc,'k')
        set(line1,'Linewidth',3)
        
        
        ylabel(strcat(caddisease,{' '},datatype))
        
        set(gca,'FontSize',GetAdjustedFontSize)
        set(gcf,'color','white')
        
        
        title(strcat(num2ordinal(index1),{' '},' Ranked Model; AICc=',num2str(AICc_best,6)))
        
        legend(strcat('Sub-epidemics=',num2str(npatches),'; C_{thr}=',num2str(onset_thr)))
        
        %legend(strcat('Num. Sub-epidemics=',num2str(npatches)))
        
        
        [npatches onset_thr]
        
        
        if (method1==0 & dist1==2)  % calculate the overdispersion factor
            
            %     [coef,ns]=getMeanVarLinear(data,totinc,6);
            %
            %     if coef>0
            %         factor1=coef;
            %     else
            
            
            % estimate dispersion in data
            binsize1=4;
            
            %[ratios,~]=getMeanVarianceRatio(smooth(data,smoothfactor1),binsize1,2);
            [ratios,~]=getMeanVarianceRatio(data,binsize1,2);
            
            index1=find(ratios(:,1)>0);
            
            factor1=mean(ratios(index1,1));
            
            factor1
            
        end
        
        
        'ABC estimates:'
        '[npatches onset_thr q]'
        
        [npatches onset_thr]
        
        factor1
        
        
        
    end
    
    
end

for j=1:cols

    subplot(rows,cols,topmodelsx(end-j+1))
    xlabel('Time');

end


save(strcat('./output/AICc_bests-modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-dateini-',datestr(datenum(caddate1),'mm-dd-yy'),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-topmodels-',num2str(topmodelsx(end)),'.mat'),'AICc_bests','-mat')


% <=====================================================================================================>
% <============================== Save file with top-ranked models' AICc metrics =======================>
% <=====================================================================================================>

performanceTop=[topmodelsx' RMSES(topmodelsx,3) relativelik_i(topmodelsx)];

T = array2table(performanceTop);
T.Properties.VariableNames(1:3) = {'i_th-ranked model','AICc','RelativeLikelihood'};
writetable(T,strcat('./output/AICc-topRanked-onsetfixed-',num2str(onset_fixed),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))



