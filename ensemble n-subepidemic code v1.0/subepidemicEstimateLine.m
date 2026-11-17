function h = subepidemicEstimateLine(estimate)
%SUBEPIDEMICESTIMATELINE Vertical dashed line at the fitted estimate.
%   h = subepidemicEstimateLine(estimate) draws a black dashed vertical line
%   at x = ESTIMATE over the full height of the current axes and returns its
%   handle. The vertical limits of the axes are kept as they were; the
%   horizontal limits widen by themselves when the estimate lies outside the
%   histogram. Call it after the histogram and with HOLD ON. Nothing is drawn
%   and h is empty when ESTIMATE is not a finite real scalar.
%
%   ESTIMATE is the value of the parameter in the ranked fit to the observed
%   data, the number shown in the title of the panel.
%
%   See also plotFit_subepidemicFramework, subepidemicParameterLabel

h = [];
if ~(isnumeric(estimate) && isreal(estimate) && isscalar(estimate) && isfinite(estimate))
    return
end
estimate = full(double(estimate));
yl = ylim;
h = plot([estimate estimate],yl,'k--');
set(h,'LineWidth',2);
ylim(yl);
end
