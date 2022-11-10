 
function [getperformance, deletetempfiles, forecastingperiod, printscreen1, weight_type1]=options_forecast


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=0; % flag or indicator variable (1/0) to calculate forecasting performance or not

deletetempfiles=1; %flag or indicator variable (1/0) to delete Forecast..mat files after use

forecastingperiod=10; %forecast horizon (number of data points ahead)

printscreen1=1;  % print plots with the results

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

weight_type1=1; % -1= equally weighted from the top models, 0=based on AICc, 1= based on relative likelihood (Akaike weights), 2=based on WISC during calibration, 3=based on WISF during forecasting performance at previous time period (week)

