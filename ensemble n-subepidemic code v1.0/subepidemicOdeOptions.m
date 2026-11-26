function options = subepidemicOdeOptions()
%SUBEPIDEMICODEOPTIONS One place for the ODE solver settings of the toolbox.
%   OPTIONS = subepidemicOdeOptions() returns the odeset structure used by
%   every fit, bootstrap refit, forecast and plot: simulateSubepidemic uses
%   it whenever no solver options are passed, and fitSubepidemicModel passes
%   it explicitly through the fit context.
%
%   Why not the odeset defaults (RelTol 1e-3, AbsTol 1e-6). Incidence is the
%   difference of two cumulative values that are 30 to 50 times larger, so a
%   relative error of 1e-3 in the cumulative curve is a much larger relative
%   error in the fitted incidence. With the defaults the AICc of the bundled
%   COVID-19 example is off by 0.5 to 2.3 units, as much as the gaps between
%   the top-ranked models, and the finite-difference gradients used by
%   fmincon and by the stationarity check are noisy.
%
%   Why tighter than the SEIR toolbox (RelTol 1e-6, AbsTol 1e-8). The
%   stationarity check of fitSubepidemicModel takes differences of the
%   objective over a step of 1e-6 of each parameter's range. In the Python
%   tests of this patch (SciPy's LSODA, not ode15s) solver noise alone gave
%   that check a value of 22 for the least-squares objective at 1e-6, above
%   its gate of 10, and 0.05 to 0.17 at 1e-8. If fitInfo.polish.run is true
%   for most of your fits, the tolerance is too loose for your data.
%
%   To change the tolerance, edit this file. Do not set Events, Mass,
%   Jacobian or Vectorized here: simulateSubepidemic owns Events and rejects
%   the others.
%
%   See also simulateSubepidemic, fitSubepidemicModel, odeset

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
end
