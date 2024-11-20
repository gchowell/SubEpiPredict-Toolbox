
function [type_GId1, mean_GI1, var_GI1] = options_Rt

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
