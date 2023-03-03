% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cumulative1,outbreakx, caddate1, cadregion, caddisease, datatype, DT, datevecfirst1, datevecend1,numstartpoints, topmodelsx, M, flag1]=options

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
global weights_obs

% <============================================================================>
% <================================ Parameters related to the data =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with the extension *.txt. The data file can contain one or more incidence curves (one per
% column in the file). Each column corresponds to the number of new cases over time for each epidemic corresponding to a different area/group.
% For instance, each column could correspond to different states in
% the U.S or countries in the world. In the options.m file, a specific data column in the file can be accessed using the parameter <outbreakx> (see below).

% if the time series file contains cumulative incidence count data, the name of the time series data file starts with "cumulative" with the
% following format:

% 'cumulative-<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt');
%  For example: 'cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt'

% Otherwise, if the time series file contains incidence data, the name of the data file follows the format:

% <cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt');
%  For example: 'daily-coronavirus-deaths-USA-05-11-2020.txt'

cumulative1=1; % flag to indicate if the data file contains cumulative incidence counts (cumulative1=1) or not (cumulative1=0)

outbreakx=52;  % identifier for the spatial area/group of interest

caddate1='05-11-2020';  % string indicating the data file date stamp in format: mm-dd-yyyy

cadregion='USA'; % string indicating the geographic region of the time series contained in the file (e.g., Georgia, USA, World, Asia, Africa, etc.)

caddisease='coronavirus'; % string indicating the name of the disease related to the time series data

datatype='deaths'; % string indicating the nature of the data (e.g., cases, deaths, hospitalizations, etc)

DT=1; % temporal resolution in days (e.g., 1=daily data, 7=weekly data, 365=yearly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end

datevecfirst1=[2020 02 27]; % 3-value date vector that specifies the date corresponding to the first data point in time series data in format: [yyy mm dd]. 

datevecend1=[2022 05 09]; % 3-value date vector that specifies the date of the most recent data file in format: [yyy mm dd].  This data file is used to assess forecast performance.

% <============================================================================>
% <============================Adjustments to data =================================>
% <============================================================================>

smoothfactor1=7; % The span of  the moving average smoothing of the case series (smoothfactor1=1 indicates no smoothing)

calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> data points where <calibrationperiod> does not exceed the length of the time series data otherwise it will use the maximum length of the data

% <=============================================================================>
% <======================= Parameter estimation and bootstrapping =========================>
% <=============================================================================>

method1=0; % Type of estimation method. See below:

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)

switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end

numstartpoints=10; % Number of initial guesses for parameter estimation procedure using MultiStart

M=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <========================= n-subepidemic growth model ===============================>
% <==============================================================================>

npatches_fixed=2; % maximum number of subepidemics considered in epidemic trajectory fit

topmodelsx=4; % Number of best fitting sub-epidemic models (based on AICc) that will be generated to derive ensemble models

if npatches_fixed==1  % if one sub-epidemic is employed, then there is only one model
    topmodelsx=1;
end

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=zeros(1,npatches_fixed)+GLM; % Sequence of subepidemic growth models considered in epidemic trajectory.

onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).

if onset_fixed==1
    if topmodelsx>npatches_fixed
        topmodelsx=npatches_fixed;
    end
end