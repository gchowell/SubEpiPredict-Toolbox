
function [outbreakx, caddate1, cadregion, caddisease, datatype, DT, datevecfirst1, datevecend1,numstartpoints, topmodelsx, M, flag1]=options

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
% The time series data file can contain one or more cumulative incidence curves (one per
% column in the file). Each column corresponds to the cumulative number of new cases for each epidemic corresponding to a different area/group.
% For instance, each column could correspond to different states in
% the U.S or countries in the world. A specific data column in the file can be accessed using the parameter <outbreakx> (see below).

% The name of the time series data file follows the following format:

% 'cumulative-<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt');
%  For example: 'cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt'

outbreakx=52;  % identifier for the spatial area of interest

caddate1='05-11-2020';  % data file time stamp in format: mm-dd-yyyy

cadregion='USA'; % string indicating the geographic region of the time series contained in the file (Georgia, USA, World, Asia, Africa, etc.)

caddisease='coronavirus'; % string indicating the name of the disease related to the time series data

datatype='deaths'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

DT=1; % temporal resolution in days (1=daily data, 7=weekly data, 365=yearly data).

datevecfirst1=[2020 02 27]; % date corresponding to the first data point in time series data in format [year_number month_number day_number]

datevecend1=[2022 05 09]; % date of the most recent data file in format [year_number month_number day_number]. This data file is accessed to assess forecast performance

% <============================================================================>
% <============================Adjustments to data =================================>
% <============================================================================>

smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series (smoothfactor1=1 indicates no smoothing)

calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> data points where <calibrationperiod> does not exceed the length of the time series data otherwise it will use the maximum length of the data

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

method1=0; % Type of estimation method. See below:

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% Pearson chi-squared=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure
%dist1=1; % error structure type (Poisson=1; NB=2)
%dist1=3; % VAR=mean+alpha*mean;
%dist1=4; % VAR=mean+alpha*mean^2;
%dist1=5; % VAR=mean+alpha*mean^d;

numstartpoints=10; % Number of initial guesses for parameter estimation procedure using MultiStart

topmodelsx=4; % Number of best fitting models (based on AICc) that will be generated to derive ensemble models

M=300; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

npatches_fixed=2; % maximum number of subepidemics considered in epidemic trajectory fit

if npatches_fixed==1  % if one sub-epidemic is employed, then there is only one model
    topmodelsx=1;
end

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=[GLM GLM]; % Sequence of subepidemic growth models considered in epidemic trajectory.

onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).

