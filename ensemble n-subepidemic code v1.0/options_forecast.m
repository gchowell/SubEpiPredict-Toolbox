function [getperformance, deletetempfiles, forecastingperiod, weight_type1] = options_forecast

% OPTIONS_FORECAST  SubEpiPredict configuration for forecasting & ensembles
%
% Overview
%   Sets the forecast horizon, whether to compute performance metrics, whether
%   to delete temporary forecast MAT files, and the weighting scheme used to
%   build an ensemble from top-ranked models.
%
% Usage
%   [getperformance, deletetempfiles, forecastingperiod, weight_type1] = options_forecast;
%
% Returns
%   getperformance     (logical) 1=compute forecast performance metrics; 0=skip
%   deletetempfiles    (logical) 1=delete temporary 'Forecast-*.mat' files; 0=keep for audit
%   forecastingperiod  (int)     Number of steps ahead to predict (horizon)
%   weight_type1       (int)     Ensemble weighting:
%                                   -1: Equal weights over the top models
%                                    0: Akaike weights, same as 1 (it used to be 1/AICc, which is not a weight)
%                                    1: Akaike weights (relative likelihoods)
%                                    2: WISC-based weights (Weighted Interval Score on calibration)
%
% Notes
%   • Ensure the weighting scheme aligns with the set of top models retained by the fit.
%
% See also
%   plotForecast_subepidemicFramework, Run_Fit_subepidemicFramework


% <==============================================================================>
% <========================== Forecasting Parameters ===========================>
% <==============================================================================>
% Define parameters for forecasting configuration.

getperformance = 1; % Flag to calculate forecasting performance metrics.
% 1 = Calculate performance metrics (e.g., error, accuracy).
% 0 = Skip performance calculations.

deletetempfiles = 1; % Flag to indicate whether temporary forecast files should be deleted.
% 1 = Delete temporary Forecast..mat files after use.
% 0 = Retain temporary Forecast..mat files for further analysis.

forecastingperiod = 4; % Forecast horizon, representing the number of time units to predict ahead.

% <==============================================================================>
% <================= Weighting Scheme for Ensemble Model =======================>
% <==============================================================================>
% Define the weighting scheme for constructing ensemble models from the top models.

weight_type1 = -1; % Type of weighting for the ensemble model:
% -1 = Equal weighting from the top models.
%  0 = Same as 1. It used to be 1/AICc, which is not a weight; the code now warns and uses Akaike weights.
%  1 = Weighted ensemble based on relative likelihood (Akaike weights).
%  2 = Weighted ensemble based on the Weighted Interval Score of the calibration period (WISC).
