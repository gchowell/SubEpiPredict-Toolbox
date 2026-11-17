function curves=AddPoissonError(yi,M,dist1,factor1,d)
%ADDPOISSONERROR Simulated incidence curves around a cumulative model curve.
%   curves = AddPoissonError(yi,M,dist1,factor1,d)
%
%   yi       cumulative model curve (a vector)
%   M        number of curves to draw
%   dist1    error structure
%              0  normal, standard deviation factor1
%              1  Poisson
%              2  negative binomial, VAR = factor1*mean
%              3  negative binomial, VAR = mean + factor1*mean
%              4  negative binomial, VAR = mean + factor1*mean^2
%              5  negative binomial, VAR = mean + factor1*mean^d
%              6  Laplace, scale factor1
%   factor1  noise scale or dispersion (not used for dist1 = 1)
%   d        power of the mean for dist1 = 5 (default 1)
%
%   curves is numel(yi)-by-M. Its first row is yi(1) in every column: the first
%   observation is not perturbed, as before. Row t+1 is drawn around the
%   increment yi(t+1) - yi(t).
%
%   The call, the distribution codes and the variance laws are those of the
%   earlier AddPoissonError.m; the body follows AddErrorStructure.m of the
%   sub-epidemic SEIR toolbox. What changed:
%
%   * A zero expected count gives 0. The earlier code computed the negative
%     binomial parameters as 0/0 and returned NaN there.
%   * Zero overdispersion (alpha = 0 for dist1 = 3, 4, 5; factor1 = 1 for
%     dist1 = 2) gives the Poisson limit. The earlier code did that only for
%     dist1 = 3 and returned NaN for the others.
%   * dist1 = 2 with factor1 < 1 is not a negative binomial. The earlier code
%     returned NaN. factor1 is an empirical variance-to-mean ratio and can fall
%     below 1, so the draw is Poisson (factor1 = 1) with warning
%     SubEpiPredict:UnderdispersedFactor. (AddErrorStructure.m stops with an
%     error there, because the SEIR toolbox floors the factor where it is
%     estimated.)
%   * Where the cumulative curve decreases, the count distributions have no
%     valid mean. The expected count is taken as 0 there. A decrease within
%     1e-6 of the largest value is ODE noise on a plateau and is passed over
%     silently; a larger one gives warning SubEpiPredict:DecreasingCurve,
%     because a model curve should not decrease. (The earlier code drew
%     around the absolute value of the increment, without a message.)
%   * Inputs are checked: a missing, empty, negative or non-finite noise scale
%     is an error that names it, not a failure further on.
%   * Draws are made one time point at a time for all M curves at once. With
%     the same random seed the numbers differ from those of the earlier code;
%     their distributions are the same.
%
%   The result never contains NaN or Inf.

if nargin<5 || isempty(d)
    d=1;
end
if nargin<4
    factor1=[];
end

if ~(isnumeric(yi) && isreal(yi) && isvector(yi) && ~isempty(yi) && all(isfinite(yi(:))))
    error('SubEpiPredict:ErrorGenerator','The cumulative curve must be a non-empty vector of finite real numbers.');
end
yi=double(full(yi(:)));
if ~(isnumeric(M) && isreal(M) && isscalar(M) && isfinite(M) && M>=0 && M==fix(M))
    error('SubEpiPredict:ErrorGenerator','The number of curves must be a nonnegative integer.');
end
if ~(isnumeric(dist1) && isreal(dist1) && isscalar(dist1) && any(dist1==0:6))
    error('SubEpiPredict:ErrorGenerator','dist1 must be an integer from 0 to 6.');
end
if dist1~=1
    if ~(isnumeric(factor1) && isreal(factor1) && isscalar(factor1) && isfinite(factor1) && factor1>=0)
        error('SubEpiPredict:ErrorGenerator', ...
            ['dist1 = %d needs factor1, the noise scale or dispersion, as one finite nonnegative ', ...
             'number; received %s. Check that method1 and dist1 in options.m match.'], ...
            dist1,localShow(factor1));
    end
    factor1=double(full(factor1));
end
if dist1==5 && ~(isnumeric(d) && isreal(d) && isscalar(d) && isfinite(d) && d>0)
    error('SubEpiPredict:ErrorGenerator','dist1 = 5 needs a finite positive power d.');
end
if dist1==2 && factor1<1
    subepidemicWarning('SubEpiPredict:UnderdispersedFactor', ...
        ['dist1 = 2 with a variance-to-mean factor of %g (below 1) is not a negative binomial. ', ...
         'Poisson noise (factor 1) is used instead.'],factor1);
    factor1=1;
end

mu=diff(yi);
if any(dist1==1:5)
    % Count means cannot be negative: a decreasing stretch has expected count 0.
    % ODE noise on a flat curve is passed over; a real decrease is reported.
    tolerance=1e-6*max(1,max(abs(yi)));
    if any(mu<-tolerance)
        subepidemicWarning('SubEpiPredict:DecreasingCurve', ...
            ['The cumulative model curve decreases at %d of %d steps (largest drop %g). ', ...
             'The expected count is taken as 0 there.'],sum(mu<-tolerance),numel(mu),max(-mu));
    end
    mu=max(mu,0);
end

curves=zeros(numel(yi),M);
curves(1,:)=yi(1);
if M==0 || isempty(mu)
    return
end

switch dist1

    case 0
        % Normal distribution
        curves(2:end,:)=repmat(mu,1,M)+factor1*randn(numel(mu),M);

    case 6
        % Laplace distribution around abs(increment), by its inverse CDF
        U=min(max(rand(numel(mu),M),eps),1-eps)-0.5;
        curves(2:end,:)=repmat(abs(mu),1,M)-factor1*sign(U).*log1p(-2*abs(U));

    otherwise
        % Poisson and negative binomial: VAR = mean*(1 + extraRatio)
        for t=1:numel(mu)
            mean1=mu(t);
            if mean1==0
                continue % a zero expected count is exactly 0
            end
            switch dist1
                case 1
                    extraRatio=0;
                case 2
                    extraRatio=factor1-1;        % VAR = factor1*mean
                case 3
                    extraRatio=factor1;          % VAR = mean + alpha*mean
                case 4
                    extraRatio=factor1*mean1;    % VAR = mean + alpha*mean^2
                case 5
                    extraRatio=factor1*mean1^(d-1); % VAR = mean + alpha*mean^d
            end
            if ~isfinite(extraRatio) || extraRatio<0
                error('SubEpiPredict:ErrorGenerator','The variance of the error structure is invalid or overflowed.');
            end
            p1=1/(1+extraRatio);
            if extraRatio==0 || p1==1
                curves(t+1,:)=poissrnd(mean1,1,M);
            else
                r1=mean1/extraRatio;
                if ~isfinite(r1) || r1<=0 || p1<=0
                    error('SubEpiPredict:ErrorGenerator','The negative binomial parameters are outside their numerical range.');
                end
                curves(t+1,:)=nbinrnd(r1,p1,1,M);
            end
        end

end

if any(~isfinite(curves(:))) || ~isreal(curves)
    error('SubEpiPredict:ErrorGenerator','The error generator produced a value that is not a finite real number.');
end
end

% -----------------------------------------------------------------------------
function s=localShow(v)
if isnumeric(v) && isempty(v)
    s='an empty value';
elseif isnumeric(v) && numel(v)<=4
    s=mat2str(double(full(v)),6);
else
    s=['a ' class(v)];
end
end
