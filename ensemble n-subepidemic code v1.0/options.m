
function [outbreakx, caddate1, cadregion, caddisease, datatype, DT, datafilename1, datevecfirst1, datevecend1,numstartpoints, topmodelsx, M, flag1]=options

% last uppdated: 11/03/22
% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 %Parameter estimation method - LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3,MLE (Neg Binomial)=4, MLE (Neg Binomial)=5

global npatches_fixed

global onset_fixed

global dist1 
global factor1
global smoothfactor1

global calibrationperiod1

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>

outbreakx=52;  % identifier for spatial area

caddate1='05-11-20';  % data file time stamp

cadregion='USA'; % string indicating the region of the time series (USA, Chile, Mexico, Nepal, etc)

caddisease='coronavirus'; % string indicating the name of the disease

datatype='deaths'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

DT=1; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
end

% Name of the file containing the cumulative time series data (rows=time, cols=regions)
datafilename1=strcat('cumulative-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-',caddate1,'.txt'); %data file with all time series across areas/regions

datevecfirst1=[2020 02 27]; % date corresponding to the first data point in time series data

datevecend1=[2022 05 09]; % date of the most recent data file which is accessed to assess forecast performance.

% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series

calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)


% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

method1=0; % Type of estimation method: 0 = LSQ

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squared=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure:

%dist1=0; % Normnal distribution to model error structure
%dist1=1; % error structure type (Poisson=1; NB=2)
%dist1=3; % VAR=mean+alpha*mean;
%dist1=4; % VAR=mean+alpha*mean^2;
%dist1=5; % VAR=mean+alpha*mean^d;

numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

topmodelsx=4; % number of best fitting models (based on AICc) that will be generated to derive ensemble models

M=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <========================= Growth model =======================================>
% <==============================================================================>

npatches_fixed=2; % maximum number of subepidemics considered in epidemic trajectory fit

if npatches_fixed==1
    topmodelsx=1;
end

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=[GLM GLM GLM]; % Sequence of subepidemic growth models considered in epidemic trajectory

onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).

