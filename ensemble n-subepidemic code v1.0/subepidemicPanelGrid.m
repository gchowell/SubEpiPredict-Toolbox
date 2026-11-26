function [rows,cols] = subepidemicPanelGrid(n)
%SUBEPIDEMICPANELGRID Rows and columns of a figure with one panel per model.
%   [rows,cols] = subepidemicPanelGrid(n) returns a grid with at least n
%   panels.
%
%   The plotting scripts built the grid from the prime factors of n:
%     one factor (n prime or 1)   n rows, 1 column
%     two factors                 factor 1 rows, factor 2 columns
%     three factors               factor 1 * factor 2 rows, factor 3 columns
%     four or more factors        factor 1 rows, factor 2 columns
%   The last case gives too few panels (2 by 2 for n = 16), so the script
%   stopped when it reached the fifth model. This function returns the same
%   grid as before whenever that grid holds n panels, which is every n with
%   at most three prime factors (1 to 15, 17 to 23, ...). Otherwise it returns
%   the squarest grid that holds n: 4 by 4 for 16, 5 by 5 for 24, 6 by 6
%   for 32.
%
%   See also plotForecast_subepidemicFramework, plotRankings_subepidemicFramework

if ~(isnumeric(n) && isreal(n) && isscalar(n) && isfinite(n) && n>=1 && n==fix(n))
    error('SubEpiPredict:PanelGrid','The number of panels must be a positive integer.');
end
n = double(n);
factors = factor(n);
if length(factors)==1
    rows = factors;
    cols = 1;
elseif length(factors)==3
    rows = factors(1)*factors(2);
    cols = factors(3);
else
    rows = factors(1);
    cols = factors(2);
end
if rows*cols<n
    cols = ceil(sqrt(n));
    rows = ceil(n/cols);
end
end
