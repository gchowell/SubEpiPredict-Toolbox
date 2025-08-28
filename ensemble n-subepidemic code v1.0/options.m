% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cumulative1, outbreakx, caddate1, cadregion, caddisease, datatype, DT, datevecfirst1, datevecend1, numstartpoints, topmodelsx, B, flag1] = options

% OPTIONS  SubEpiPredict configuration for fitting n-sub-epidemic models
%
% Overview
%   Defines data tags, temporal resolution, smoothing/calibration settings,
%   estimation/error models, and the search space for n-sub-epidemic growth
%   trajectories. Used by the fitting scripts to calibrate candidate models
%   (1..n patches) and rank them (e.g., by AICc).
%
% Usage
%   [cumulative1, outbreakx, caddate1, cadregion, caddisease, datatype, ...
%    DT, datevecfirst1, datevecend1, numstartpoints, topmodelsx, B, flag1] = options;
%
% Returns
%   cumulative1      (logical) 1 if input file stores cumulative counts; 0 if incidence
%   outbreakx        (int)     COLUMN INDEX of the focal area/group in the input matrix
%   caddate1         (char)    Datestamp in filenames, format 'mm-dd-yyyy'
%   cadregion        (char)    Geographic tag (e.g., 'USA')
%   caddisease       (char)    Disease tag (e.g., 'coronavirus')
%   datatype         (char)    'cases' | 'deaths' | 'hospitalizations' | …
%   DT               (int)     Time step: 1=daily, 7=weekly, 365=yearly
%   datevecfirst1    (1×3 int) First date present [yyyy mm dd]
%   datevecend1      (1×3 int) Last date present  [yyyy mm dd]
%   numstartpoints   (int)     MultiStart initial points for optimization
%   topmodelsx       (int)     # top-ranked models to keep (by AICc)
%   B                (int)     Bootstrap replicates for parameter/forecast uncertainty
%   flag1            (int)     Growth-kernel code used in the trajectory (see below)
%
% Globals used/set
%   method1            Estimator:
%                        0=LSQ, 1=MLE Poisson, 2=Pearson χ², 3/4/5=MLE NegBin (see below)
%   npatches_fixed     Max number of sub-epidemic components (patches)
%   onset_fixed        1=fixed onset times across patches; 0=asynchronous allowed
%   dist1              Observation/error model code (mapped below)
%   smoothfactor1      Moving-average span for pre-smoothing (1=no smoothing)
%   calibrationperiod1 # of most recent points used for calibration
%
% Input data & filename patterns (in ./input)
%   If cumulative1==1 → file name MUST start with 'cumulative-':
%     'cumulative-<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt'
%   Else (incidence):
%     '<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt'
%   where <cadtemporal> ∈ {'daily','weekly','yearly'}; columns = areas/groups.
%
% Estimation & error models
%   method1 chooses the estimator; dist1 encodes/weights the observation model:
%     method1 = 0  LSQ
%         choose dist1 ∈ {0,1,2}:
%           0: Normal (homoscedastic LSQ)
%           1: Poisson-like weighting (var≈mean; LSQ variant)
%           2: NegBin-like weighting var=factor1·mean (factor1 estimated empirically)
%     method1 = 1  MLE Poisson                  → dist1 := 1 (automatic)
%     method1 = 2  Pearson χ² (GOF diagnostic)  → dist1 unchanged (no auto-mapping)
%     method1 = 3  MLE NegBin var = mean + α·mean    → dist1 := 3 (automatic)
%     method1 = 4  MLE NegBin var = mean + α·mean^2  → dist1 := 4 (automatic)
%     method1 = 5  MLE NegBin var = mean + α·mean^d  → dist1 := 5 (automatic)
%
% n-sub-epidemic growth kernels (flag1 codes)
%   GGM=0 (Generalized Growth), GLM=1 (Logistic), GRM=2 (Generalized Richards),
%   LM=3 (Linear), RICH=4 (Richards). Set flag1 to one of the above.
%
% Model search & constraints
%   npatches_fixed caps the number of sub-epidemic components considered.
%   If onset_fixed==1 and topmodelsx > npatches_fixed, topmodelsx is reduced to npatches_fixed.
%   If npatches_fixed==1, topmodelsx is forced to 1.
%
% Notes
%   • Smoothing via smoothfactor1 helps stabilize noisy incidence series.
%   • Calibration uses the last calibrationperiod1 points of the series.
%
% See also
%   Run_Fit_subepidemicFramework, plotFit_subepidemicFramework, plotRankings_subepidemicFramework


% <============================================================================>
% <=================== Declare Global Variables ===============================>
% <============================================================================>
% Declare global variables used for parameter estimation and model fitting.
global method1          % Parameter estimation method: LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3, 4, 5
global npatches_fixed   % Maximum number of subepidemics considered in the trajectory
global onset_fixed      % Flag to indicate if the onset timing of subepidemics is fixed
global dist1            % Type of error structure
global smoothfactor1    % Smoothing factor for time series data
global calibrationperiod1 % Calibration period for the model

% <============================================================================>
% <========================= Dataset Properties ===============================>
% <============================================================================>
% The input folder contains a time series data file in *.txt format. The file can contain:
% - One or more incidence curves (columns represent spatial areas/groups).
% - Columns contain the number of new cases over time for each group/region.

% If the data file contains cumulative incidence counts, its name starts with "cumulative":
% Format: 'cumulative-<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt'
% Example: 'cumulative-daily-coronavirus-deaths-USA-05-11-2020.txt'

% Otherwise, for incidence data:
% Format: '<cadtemporal>-<caddisease>-<datatype>-<cadregion>-<caddate1>.txt'
% Example: 'daily-coronavirus-deaths-USA-05-11-2020.txt'

cumulative1 = 1;        % Flag: 1 if data file contains cumulative counts; 0 otherwise
outbreakx = 52;         % Identifier for the spatial area/group of interest
caddate1 = '05-11-2020'; % Data file date stamp (format: mm-dd-yyyy)
cadregion = 'USA';      % Geographic region of the time series (e.g., Georgia, USA, Asia)
caddisease = 'coronavirus'; % Name of the disease (e.g., coronavirus)
datatype = 'deaths';    % Type of data (e.g., cases, deaths, hospitalizations)
DT = 1;                 % Temporal resolution (1=daily, 7=weekly, 365=yearly)

% Temporal resolution description
if DT == 1
    cadtemporal = 'daily';
elseif DT == 7
    cadtemporal = 'weekly';
elseif DT == 365
    cadtemporal = 'yearly';
end

% Dates corresponding to the time series
datevecfirst1 = [2020 02 27]; % Date of the first data point [yyyy mm dd]
datevecend1 = [2022 05 09];   % Most recent data file date [yyyy mm dd]

% <============================================================================>
% <============================ Data Adjustments ==============================>
% <============================================================================>
% Smoothing and calibration settings for time series data.

smoothfactor1 = 7;      % Moving average smoothing span (1=no smoothing)
calibrationperiod1 = 90; % Number of most recent data points used for calibration
% If calibration period exceeds the time series length, the maximum length is used.

% <============================================================================>
% <================== Parameter Estimation and Bootstrapping ==================>
% <============================================================================>
% Settings for parameter estimation and error structure.

method1 = 0;            % Parameter estimation method:
% 0 = Nonlinear least squares (LSQ)
% 1 = MLE Poisson
% 3 = MLE Negative Binomial (VAR=mean+alpha*mean)
% 4 = MLE Negative Binomial (VAR=mean+alpha*mean^2)
% 5 = MLE Negative Binomial (VAR=mean+alpha*mean^d)

dist1 = 0;              % Error structure type:
% 0 = Normal distribution
% 1 = Poisson error structure
% 2 = Negative Binomial (VAR=factor1*mean)
% 3 = MLE Negative Binomial (VAR=mean+alpha*mean)
% 4 = MLE Negative Binomial (VAR=mean+alpha*mean^2)
% 5 = MLE Negative Binomial (VAR=mean+alpha*mean^d)

% Automatically set dist1 based on method1
switch method1
    case 1
        dist1 = 1;
    case 3
        dist1 = 3;
    case 4
        dist1 = 4;
    case 5
        dist1 = 5;
end

numstartpoints = 10;    % Number of initial guesses for parameter estimation (MultiStart)
B = 300;                % Number of bootstrap realizations for uncertainty characterization

% <============================================================================>
% <================= n-Subepidemic Growth Model Settings ======================>
% <============================================================================>
% Configuration for subepidemic growth models.

npatches_fixed = 2;     % Maximum number of subepidemics in the model fit
topmodelsx = 4;         % Number of best-fitting subepidemic models (based on AICc)

% Adjust number of top models if a single subepidemic is used
if npatches_fixed == 1
    topmodelsx = 1;
end

% Growth model types
GGM = 0;  % Generalized Growth Model
GLM = 1;  % Logistic Model
GRM = 2;  % Generalized Richards Model
LM = 3;   % Linear Model
RICH = 4; % Richards Model

flag1 = GLM;            % Growth model sequence used in trajectory fitting
onset_fixed = 0;        % Fix onset timing of subepidemics (1=fixed, 0=not fixed)

% Ensure number of top models does not exceed subepidemics when onset is fixed
if onset_fixed == 1 && topmodelsx > npatches_fixed
    topmodelsx = npatches_fixed;
end
