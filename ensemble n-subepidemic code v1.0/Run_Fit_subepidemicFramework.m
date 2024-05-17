% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function Run_Fit_subepidemicFramework(outbreakx_pass,caddate1_pass)

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

clear global

global method1 %Parameter estimation method - LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3,MLE (Neg Binomial)=4, MLE (Neg Binomial)=5
global npatches_fixed
global onset_fixed
global dist1
global factor1
global smoothfactor1
global calibrationperiod1

% <============================================================================>
% <================== Load the parameter values ===============================>
% <============================================================================>

[cumulative1_INP, outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options


% <============================================================================>
% <================================ Datasets properties =======================>
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

caddisease=caddisease_INP; % string indicating the name of the disease

datatype=datatype_INP; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data, 365= yearly data).

cumulative1=cumulative1_INP;

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end

% Name of the file containing the cumulative time series data (rows=time, cols=regions)

if cumulative1==1
    datafilename1=strcat('cumulative-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-',caddate1,'.txt'); %data file with all time series across areas/regions
else
    datafilename1=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-',caddate1,'.txt'); %data file with all time series across areas/regions
end

% Name of the file for the adjusted incidence data file for a specific region and
% after removing early zeros.

datafilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1,'.txt'); % specific adjusted data file for a given region/area

datevecfirst1=datevecfirst1_INP; % date corresponding to the first data point in time series data

datevecend1=datevecend1_INP;

% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

%smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)


% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squared=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

% % Define dist1 which is the type of error structure:
% switch method1
%     
%     case 0
%         
%         dist1=0; % Normnal distribution to model error structure
%        
%         %dist1=1; % error structure type (Poisson=1; NB=2)
%         
%         %factor1=1; % scaling factor for VAR=factor1*mean
%         
%         
%     case 3
%         dist1=3; % VAR=mean+alpha*mean;
%         
%     case 4
%         dist1=4; % VAR=mean+alpha*mean^2;
%         
%     case 5
%         dist1=5; % VAR=mean+alpha*mean^d;
%         
% end

numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

topmodelsx=topmodelsx_INP; % number of best fitting models (based on AICc) that will be generated to derive ensemble models

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <========================= Mathematical model ======================================>
% <==============================================================================>

%npatches_fixed=2; % maximum number of subepidemics considered in epidemic trajectory fit

if npatches_fixed==1
    topmodelsx=1;
end

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory

%flag1=[RICH RICH];

%onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).

% <===========================================================================================================>
% <====== Check that the number of estimated parameters is smaller than the number of data points= ===========>
% <===========================================================================================================>

numparams=get_nparams(method1,dist1,npatches_fixed,flag1,1,onset_fixed);

numparams
calibrationperiod1

if numparams>=calibrationperiod1

    error("Number of estimated parameters should be smaller than the calibration period. Consider reducing <npatches_fixed> or increasing the length of the calibration period.")

end


% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>

%pause(10*rand)
data=load(strcat('./input/',datafilename1)); % load time series dataset (rows=time, cols=regions/areas)

if isempty(data)
    error('The dataset is empty')
end

dataprov=data';

timelags=0; % keeps track of the epidemic onset timing (nonzero case)

for outbreak1=outbreakx
    
    close all

    data1=dataprov(outbreak1,:)';

    if cumulative1==1

        data1=[data1(1);diff(data1)]; % Incidence curve

    end

    clear dataprov
    clear data
    
    % <==================================================================================>
    % <============ select length of the calibration period =============================>
    % <==================================================================================>
    
    % select length of the calibration period
    if calibrationperiod1<length(data1)
        
        data1(1:1:end-calibrationperiod1)=0;
        
    end
    
    %plot(data1)
    %hold on
    
    % <=============================================================================================>
    % <========== remove early zeros until a non-zero data point is found in the time series===============================>
    % <=============================================================================================>
    
    index1=find(data1>0);
    
    data1=data1(index1(1):end,1);
    
    timelags=index1(1)-1;
    


    %     if 0
    %         % remove early low testing period and start from the monotonic increasing trend
    %
    %         [max2, index2]=max(data1);
    %
    %         [min1, index1]=min(data1(1:index2));
    %         [max2, index2]=max(data1(index1:end));
    %
    %         index2=index1+index2-1;
    %
    %         [index1, min1]=find(data1(1:index2)==min1);
    %         index1=index1(length(index1));
    %
    %
    %         if index1<index2
    %             data1=data1(index1:end);
    %
    %             timelags=timelags+index1(1)-1;
    %
    %         else
    %             data1=[];
    %             length(data1)
    %             outbreak1
    %             'small number of data points..1'
    %
    %             continue
    %         end
    %
    %         if data1(1)==0
    %             data1=data1(2:end,1);
    %             timelags=timelags+1;
    %         end
    %     end
    
    outbreakx
    
    caddate1
    
    data1
    
    % <=============================================================================================>
    % <================ If time series is too short, do not proceed further ==========================================>
    % <=============================================================================================>
    
    if length(data1)<10
        length(data1)
        'small number of data points..'
        continue
    end
    
    
    if (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure
        
        binsize1=7;
        
        [ratios,~]=getMeanVarianceRatio(data1,binsize1,2);
        
        index1=find(ratios(:,1)>0);
        
        factor1=mean(ratios(index1,1));
        
        factor1
        
    end
    
    
    % <==============================================================================================>
    % <================ Save adjusted incidence data file for further analysis ========================================>
    % <==============================================================================================>
    
    data1=[(0:1:length(data1)-1)' data1];
    
    %save(strcat('./output/',datafilename2),'data1','-ascii')
    
    data=data1;
    
    t_window=1:length(data(:,1));
    
    data1=data(t_window',:)
    
    plot(data1(:,1),data1(:,2),'ko')
    
    % <=======================================================================================>
    % <===================== Derive the best fitting sub-epidemic models ==================================>
    % <=======================================================================================>
    
    [RMSESx,PS,npatches,onset_thr,~]=fittingModifiedLogisticFunctionPatchABC(datafilename2,data1,DT,t_window,M,flag1,numstartpoints);

    AICmin=RMSESx(1,3);

    relativelik_i=exp((AICmin-RMSESx(:,3))/2);
    
    
    if 1  %plot relative likelihoods of the models
        
        figure(99)
        
        subplot(1,3,1)
        line1=plot(RMSESx(1:topmodelsx,3),'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('AICc')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        subplot(1,3,2)
        
        line1=plot(relativelik_i(1:topmodelsx),'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('Relative likelihood')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        subplot(1,3,3)
        
        %deltas=AICc_bests-AICc_bests(1);
        
        evidenceRatio=1./relativelik_i(1:topmodelsx);
        
        line1=plot(evidenceRatio,'ko-')
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('Evidence ratio')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
    end
    
    index1=find(evidenceRatio<100); %check top models with evidence ratio < 100
    
    %if length(index1)>2 %if there are more than two top models with evidence ratio < 100
        
    %    topmodelsx=length(index1);
        
    %else
        
    %    topmodelsx=2;
        
    %end
    
    
    % <===========================================================================================>
    % <=============== Derive uncertainty for the <topmodelsx> best fitting models and save results =====================>
    % <===========================================================================================>
      
   
    for rank1=1:topmodelsx
        
        [Phatss,npatches,onset_thr,curves,bestfit,data1,P0,AICc_best,factor1,d]=fittingModifiedLogisticFunctionPatchMultiple(RMSESx,PS,data1,DT,t_window,M,flag1,numstartpoints,rank1);
        
        cadfilename1=strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename2(1:end-4),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat');
        
        save(cadfilename1,'-mat')
        
        strcat('rank-',num2str(rank1),'-complete!')
    end
    
end

clear RMSESx PS P
clear Phatss curves bestfit data1 P0 AICc_best d

