% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function Run_Fit_subepidemicFramework(outbreakx_pass,caddate1_pass)

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

% Reset the toolbox's own global variables only. A bare "clear global" also
% deleted every other global variable in the MATLAB session.
clear global method1 npatches_fixed onset_fixed dist1 factor1 smoothfactor1 calibrationperiod1
clear global invasions timeinvasions Cinvasions npatches onset_thr flag1 I0
clear global timevect ydata yfit data LBe UBe

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

[cumulative1_INP, outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options;


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


% <============================================================================>
% <== Check the settings before anything else: every problem is reported at once ==>
% <============================================================================>

checkedSettings=struct();
checkedSettings.method1=method1;
checkedSettings.dist1=dist1;
checkedSettings.flag1=flag1_INP;
checkedSettings.npatches_fixed=npatches_fixed;
checkedSettings.onset_fixed=onset_fixed;
checkedSettings.smoothfactor1=smoothfactor1;
checkedSettings.calibrationperiod1=calibrationperiod1;
checkedSettings.numstartpoints=numstartpoints_INP;
checkedSettings.topmodelsx=topmodelsx_INP;
checkedSettings.B=M_INP;
checkedSettings.DT=DT_INP;
checkedSettings.cumulative1=cumulative1_INP;
checkedSettings.outbreakx=outbreakx;
subepidemicCheckSettings(checkedSettings);

% Revision i: the settings in one line (the 13 outputs of options.m are no
% longer displayed one by one), and a fresh list of warnings for this run.
subepidemicWarningLog('start');
fprintf(['Settings: method1 = %d, dist1 = %d, flag1 = %d, up to %d sub-epidemic(s), onset_fixed = %d, ', ...
    'smoothfactor1 = %d, calibrationperiod1 = %d, %d start points, %d top models, B = %d.\n'], ...
    method1,dist1,flag1_INP,npatches_fixed,onset_fixed,smoothfactor1,calibrationperiod1, ...
    numstartpoints_INP,topmodelsx_INP,M_INP);
clear checkedSettings

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

numparams;
calibrationperiod1;

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

% Several areas: each area is fitted and saved as if the driver had been
% called for that area alone. Before, (1) the data matrix was cleared after
% the first area, so the second one failed; (2) the output file names were
% built once from the whole vector outbreakx, so every area would have written
% to the same files; (3) a column vector was not looped over at all.
outbreakList=outbreakx(:)';
areasSkipped=[]; % areas with too few data points (revision h)
topmodelsRequested=topmodelsx; % the setting; topmodelsx can be lowered for one area (revision h)

% Every requested area must be a column of the data file; checked for all
% areas now, not after the earlier ones have been fitted.
if any(outbreakList>size(dataprov,1))
    error('SubEpiPredict:InvalidSettings', ...
        'outbreakx = %s, but the data file %s has %d column(s).',mat2str(outbreakList),datafilename1,size(dataprov,1));
end

for outbreak1=outbreakList

    close all

    outbreakx=outbreak1; % the saved workspace and the file names refer to this area

    topmodelsx=topmodelsRequested;

    % Random seed (setting rngSeed) and generator state, per area, so that an
    % area fitted in a call with several areas equals the same area fitted
    % alone. rngInfo is saved with the results of each ranked model.
    rngInfo=subepidemicSeed(sprintf('Run_Fit_subepidemicFramework, area %d',outbreak1));

    datafilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreak1),'-',caddate1,'.txt');

    data=load(strcat('./input/',datafilename1));
    dataprov=data';

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
    
    outbreakx;
    
    caddate1;
    
    data1;
    
    % <=============================================================================================>
    % <================ If time series is too short, do not proceed further ==========================================>
    % <=============================================================================================>
    
    if length(data1)<10 && npatches_fixed>1

        % Revision h: say so. Before, the area was skipped with a printed
        % string only and the driver ended as if it had been fitted.
        subepidemicWarning('SubEpiPredict:AreaSkipped', ...
            ['Area %d has %d data points after the calibration window and the leading zeros are removed. ', ...
             'At least 10 are needed with npatches_fixed > 1. The area was skipped: nothing was fitted or saved for it.'], ...
            outbreak1,length(data1));
        areasSkipped(end+1)=outbreak1; %#ok<AGROW>
        length(data1);
        'small number of data points..';
        continue
    end
    
    
    if (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure
        
        binsize1=7;
        
        [ratios,~]=getMeanVarianceRatio(data1,binsize1,2);
        
        index1=find(ratios(:,1)>0);
        
        factor1=mean(ratios(index1,1));
        
        factor1;
        
    end
    
    
    % <==============================================================================================>
    % <================ Save adjusted incidence data file for further analysis ========================================>
    % <==============================================================================================>
    
    data1=[(0:1:length(data1)-1)' data1];
    
    %save(strcat('./output/',datafilename2),'data1','-ascii')
    
    data=data1;
    
    t_window=1:length(data(:,1));
    
    data1=data(t_window',:);
    
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    nexttile(1)
    plot(data1(:,1),data1(:,2),'ko')
    
    % <=======================================================================================>
    % <===================== Derive the best fitting sub-epidemic models ==================================>
    % <=======================================================================================>
    
    [RMSESx,PS,npatches,onset_thr,~]=fittingModifiedLogisticFunctionPatchABC(datafilename2,data1,DT,t_window,M,flag1,numstartpoints);

    % Revision h. The screen can retain fewer models than topmodelsx (few
    % onset thresholds, candidates too large for the series, rejected fits).
    % The lines below indexed 1:topmodelsx and stopped with an index error
    % after the whole screen.
    if topmodelsx>size(RMSESx,1)
        subepidemicWarning('SubEpiPredict:TopModelsCapped', ...
            ['options.m asks for %d top-ranked models, but the screen of area %d retained %d. ', ...
             'Only those %d are bootstrapped.'],topmodelsx,outbreak1,size(RMSESx,1),size(RMSESx,1));
        topmodelsx=size(RMSESx,1);
    end

    AICmin=RMSESx(1,3);

    relativelik_i=exp((AICmin-RMSESx(:,3))/2);

    % Revision l: evidence ratios of the top models against ranked model 1,
    % with a warning for those the data hardly support (report only; the
    % number of top models does not change). See subepidemicEvidenceRatios.m.
    fitSettingsArea=subepidemicFitSettings();
    % Revision m1: what made this fit (revision and settings); saved with every
    % ranked file and compared by subepidemicRollingForecast.m before a fit is reused.
    fitSignature=subepidemicFitSignature(); %#ok<NASGU>
    evidenceRatio=subepidemicEvidenceRatios(RMSESx(:,3),topmodelsx,outbreak1,fitSettingsArea.evidenceRatioLimit);
    fprintf('Evidence ratios of ranked models 1 to %d against ranked model 1: %s\n',topmodelsx, ...
        strtrim(sprintf('%.3g  ',full(evidenceRatio))));
    
    
    if 1  %plot relative likelihoods of the models
        
        figure(99)
        
        tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

        nexttile(1)
        line1=plot(RMSESx(1:topmodelsx,3),'ko-');
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('AICc')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        nexttile(2)
        
        line1=plot(relativelik_i(1:topmodelsx),'ko-');
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('Relative likelihood')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
        nexttile(3)
        
        %deltas=AICc_bests-AICc_bests(1);
        
        % evidenceRatio is computed above (subepidemicEvidenceRatios.m)
        
        line1=plot(evidenceRatio,'ko-');
        
        set(line1,'LineWidth',2)
        
        xlabel('i_{th} ranked model')
        ylabel('Evidence ratio')
        set(gca,'FontSize',GetAdjustedFontSize);
        set(gcf,'color','white')
        
        
    end
    
    % Before revision l a line here found the top models with an evidence ratio
    % below 100 and did nothing with them (the code that would have reduced
    % topmodelsx was commented out). The ratio is now reported above.


    % <===========================================================================================>
    % <=============== Derive uncertainty for the <topmodelsx> best fitting models and save results =====================>
    % <===========================================================================================>
      
   
    for rank1=1:topmodelsx
        
        [Phatss,npatches,onset_thr,curves,bestfit,data1,P0,AICc_best,factor1,d,bootstrapDiagnostics]=fittingModifiedLogisticFunctionPatchMultiple(RMSESx,PS,data1,DT,t_window,M,flag1,numstartpoints,rank1);
        
        cadfilename1=strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename2(1:end-4),'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat');
        
        save(cadfilename1,'-mat')
        
        fprintf('Ranked model %d complete.\n',rank1);
    end
    
end

if ~isempty(areasSkipped)
    fprintf('Areas skipped for too few data points (nothing was saved for them): %s\n',mat2str(areasSkipped));
end

clear RMSESx PS P
clear Phatss curves bestfit data1 P0 AICc_best d

% Revision i: list the toolbox's warnings of this run once more (subepidemicWarningLog.m).
subepidemicWarningLog('summary','Run_Fit_subepidemicFramework');
