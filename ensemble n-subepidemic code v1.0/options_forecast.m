 
function [getperformance, deletetempfiles, forecastingperiod, printscreen1, weight_type1]=options_forecast


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance metrics or not

deletetempfiles=1; %flag or indicator variable (1/0) to indicate whether we wan to delete Forecast..mat files after use

forecastingperiod=30; % forecast horizon (number of time units ahead)

printscreen1=1;  % flag or indicator variable (1/0) to indicate whether we weant to print plots with the results

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

weight_type1=1; % -1= equally weighted from the top models, 0= weighted ensemble based on AICc, 1= weighted ensemble based on relative likelihood (Akaike weights), 
% 2=weighted ensemble based on WISC during calibration, 3=weighted ensemble based on WISF during forecasting performance at previous time period (week)

