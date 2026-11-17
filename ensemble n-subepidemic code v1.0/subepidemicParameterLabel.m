function label = subepidemicParameterLabel(name,estimate,lower,upper)
%SUBEPIDEMICPARAMETERLABEL Title text for one parameter: estimate and interval.
%   label = subepidemicParameterLabel(name,estimate,lower,upper) returns
%       '<name>=<estimate>(95%CI:<lower>,<upper>)'
%   for example 'r_2=24.2(95%CI:7.63,102)' or 'K_1=45600(95%CI:39600,52400)'.
%
%   ESTIMATE is the value of the parameter in the ranked fit (the fit to the
%   observed data). LOWER and UPPER are the 2.5% and 97.5% quantiles of the
%   bootstrap refits. The mean of the refits is not used: for a skewed
%   parameter it differs from the estimate.
%
%   Each number is rounded to three significant digits and written without
%   an exponent when its magnitude is between 0.001 and 10 million:
%   0.348, 24.2, 102, 45600. Trailing zeros are dropped (0.46, 1).
%
%   See also plotFit_subepidemicFramework

if ~(ischar(name) && (isrow(name) || isempty(name)))
    error('SubEpiPredict:InvalidLabel','name must be a character row.');
end
values = {estimate,lower,upper};
for k = 1:3
    v = values{k};
    if ~(isnumeric(v) && isreal(v) && isscalar(v))
        error('SubEpiPredict:InvalidLabel','estimate, lower and upper must be real scalars.');
    end
end
label = sprintf('%s=%s(95%%CI:%s,%s)',name,localNumber(estimate),localNumber(lower),localNumber(upper));
end

% -----------------------------------------------------------------------------
function s = localNumber(x)
% Three significant digits; plain notation for 0.001 <= |x| < 1e7.
x = full(double(x));
if ~isfinite(x)
    s = sprintf('%g',x);
    return
end
if x==0
    s = '0';
    return
end
e = floor(log10(abs(x)));
y = round(x/10^(e-2))*10^(e-2);          % three significant digits
if abs(y)>=1e7 || abs(y)<1e-3
    s = sprintf('%.3g',y);
elseif abs(y)>=1000
    s = sprintf('%d',round(y));
else
    s = sprintf('%.3g',y);
end
end
