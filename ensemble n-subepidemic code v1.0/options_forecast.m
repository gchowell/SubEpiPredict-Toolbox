function [getperformance, deletetempfiles, forecastingperiod, weight_type1] = options_forecast

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

forecastingperiod = 30; % Forecast horizon, representing the number of time units to predict ahead.

% <==============================================================================>
% <================= Weighting Scheme for Ensemble Model =======================>
% <==============================================================================>
% Define the weighting scheme for constructing ensemble models from the top models.

weight_type1 = 1; % Type of weighting for the ensemble model:
% -1 = Equal weighting from the top models.
%  0 = Weighted ensemble based on Akaike Information Criterion corrected (AICc).
%  1 = Weighted ensemble based on relative likelihood (Akaike weights).
%  2 = Weighted ensemble based on the Weighted Interval Score of the calibration period (WISC).
