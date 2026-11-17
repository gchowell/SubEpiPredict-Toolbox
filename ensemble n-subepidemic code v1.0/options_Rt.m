
function [type_GId1, mean_GI1, var_GI1, minGIcoverage] = options_Rt

% OPTIONS_RT  SubEpiPredict configuration for effective reproduction number (Rt)
%
% Overview
%   Specifies the generation-interval (GI) model and parameters used to
%   compute Rt from calibrated sub-epidemic trajectories (and optionally to
%   project Rt). Keep GI units consistent with the data time step.
%
% Usage
%   [type_GId1, mean_GI1, var_GI1, minGIcoverage] = options_Rt;
%
% Returns
%   type_GId1  (int)     GI distribution family:
%                          1=Gamma, 2=Exponential, 3=Delta (fixed interval)
%   mean_GI1   (double)  Mean GI (same time units as the data, e.g., days)
%   var_GI1    (double)  Variance of the GI (units^2)
%   minGIcoverage (double) Rt is reported only where the curve so far covers at
%                        least this share of the GI distribution (0 to 1).
%                        Earlier values are blank; 0 turns the blanking off.
%
% Notes
%   • Choose GI parameters from literature appropriate to your pathogen/context.
%   • Verify consistency of time units: if DT=7 (weekly), GI parameters should be in weeks.
%
% See also
%   plotReproductionNumber, plotForecast_subepidemicFramework


% <=======================================================================================>
% <========================== Reproduction Number Parameters =============================>
% <=======================================================================================>
% This section defines parameters related to the reproduction number (Rt) calculation.

type_GId1 = 1; % Type of generation interval distribution:
               % 1 = Gamma distribution.
               % 2 = Exponential distribution.
               % 3 = Delta distribution (fixed generation interval).

mean_GI1 = 11.7; % Mean of the generation interval distribution (in time units, e.g., days).

var_GI1 = 3.7;   % Variance of the generation interval distribution (in time units squared).

minGIcoverage = 0.99; % Rt is blank (NaN) until the curve covers this share of the
                      % generation-interval distribution. Before that the renewal
                      % equation misses infectors from before the first data point
                      % and overstates Rt (see get_Rt.m). 0 = no blanking.
