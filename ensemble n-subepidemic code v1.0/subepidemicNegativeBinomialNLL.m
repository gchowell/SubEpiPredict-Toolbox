function nll = subepidemicNegativeBinomialNLL(y,mu,alpha,d,method)
%SUBEPIDEMICNEGATIVEBINOMIALNLL Vectorized legacy NB fitting objective.
%   NLL = subepidemicNegativeBinomialNLL(Y,MU,ALPHA,D,METHOD) implements
%     METHOD 3: variance = mu + alpha*mu
%     METHOD 4: variance = mu + alpha*mu^2
%     METHOD 5: variance = mu + alpha*mu^d
%   It omits sum(gammaln(Y+1)), just as BOTH original objective wrappers do.
%   Row/column vectors are accepted; Y and MU must have equal lengths.
%
%   IMPORTANT: the former loop 0:(Y(i)-1) used floor(Y(i)) product factors
%   for ordinary nonnegative fractional Y, but the other terms used Y itself.
%   This performance patch explicitly preserves that floor-product convention:
%       gammaln(r + floor(Y)) - gammaln(r).
%   It does NOT silently substitute a gamma-extended fractional likelihood or
%   round the entire series. Fractional data with this convention do not define
%   a standard NB count likelihood. At floating-point integer boundaries, the
%   floor convention is explicit; platform-dependent colon endpoint rounding
%   is not emulated. See NB_LIKELIHOOD_PATCH.md.
%
%   MU must be strictly positive: the calling wrappers retain their existing
%   zero-to-0.001 replacement and all-zero-trajectory penalty. Invalid numeric
%   parameter/data values return +Inf, not NaN/complex values. Shape/method
%   errors throw descriptive exceptions. At ALPHA=0, integer Y use the exact
%   Poisson limit; fractional Y return +Inf (their legacy objective diverges).
%
%   No count-dependent loops or additional toolboxes are required. Log-domain
%   parameterization avoids powers overflowing. For large r, a rearranged
%   Stirling difference avoids cancellation in gammaln(r+n)-gammaln(r).
%   Reference: NIST DLMF 5.11.1; MATLAB documentation for gammaln/log1p/expm1.

if ~(isnumeric(y) && isvector(y) && ~isempty(y) && ...
        isnumeric(mu) && isvector(mu) && numel(mu)==numel(y))
    error('SubEpiPredict:NBInputShape', ...
        'Y and MU must be nonempty numeric vectors of equal length.');
end
if ~(isnumeric(method) && isreal(method) && isscalar(method) && ...
        any(method==[3 4 5]))
    error('SubEpiPredict:NBMethod','NB method must be 3, 4, or 5.');
end
nll = Inf;
if ~isreal(y) || ~isreal(mu) || ...
        ~(isnumeric(alpha) && isreal(alpha) && isscalar(alpha) && ...
          isfinite(alpha) && alpha>=0) || ...
        ~(isnumeric(d) && isreal(d) && isscalar(d) && isfinite(d))
    return
end
y = double(y(:));
mu = double(mu(:));
if any(~isfinite(y) | y<0 | ~isfinite(mu) | mu<=0)
    return
end
alpha = double(alpha);
d = double(d);
n = floor(y);
logMu = log(mu);

if alpha==0
    if all(y==n)
        nll = sum(mu-y.*logMu);
        if ~isfinite(nll), nll = Inf; end
    end
    return
end

logAlpha = log(alpha);
switch method
    case 3
        logBeta = repmat(logAlpha,size(mu));
        logR = logMu-logAlpha;
    case 4
        logBeta = logAlpha+logMu;
        logR = repmat(-logAlpha,size(mu));
    case 5
        logBeta = logAlpha+(d-1).*logMu;
        logR = -logAlpha+(2-d).*logMu;
end
if any(~isfinite(logBeta) | ~isfinite(logR))
    return
end

% beta = mu/r. Stable log(1+beta), log(beta/(1+beta)), and r*log(1+beta).
% Computing the last quantity as (mu/beta)*log(1+beta) directly can give
% Inf*0 or 0/0. The two branches below retain the limiting value mu.
logOnePlusBeta = max(logBeta,0)+log1p(exp(-abs(logBeta)));
logFailure = zeros(size(y));
penalty = zeros(size(y));
smallBeta = logBeta<0;
if any(smallBeta)
    b = exp(logBeta(smallBeta));
    ratio = ones(size(b));
    nonzero = b>0;
    ratio(nonzero) = log1p(b(nonzero))./b(nonzero);
    logFailure(smallBeta) = logBeta(smallBeta)-log1p(b);
    penalty(smallBeta) = mu(smallBeta).*ratio;
end
if any(~smallBeta)
    logFailure(~smallBeta) = -log1p(exp(-logBeta(~smallBeta)));
    penalty(~smallBeta) = exp(logR(~smallBeta)+ ...
        log(logOnePlusBeta(~smallBeta)));
end

logTerms = zeros(size(y));
largeR = logR>=log(16);
if any(~largeR)
    % Gamma recurrence removes gammaln(r)'s singularity as r approaches 0:
    % log Gamma(r+n)-log Gamma(r) = log(r)+log Gamma(r+n)-log Gamma(r+1).
    % It remains usable when exp(logR) underflows; n=0 has an empty product.
    ii = find(~largeR);
    ni = n(ii);
    r = exp(logR(ii));
    rising = zeros(size(ni));
    positive = ni>0;
    rising(positive) = logR(ii(positive))+ ...
        gammaln(r(positive)+ni(positive))-gammaln(1+r(positive));
    logTerms(ii) = rising+y(ii).*logFailure(ii)-penalty(ii);
end
if any(largeR)
    ii = find(largeR);
    ni = n(ii);
    correction = scaledLogRisingFactorial(ni,logR(ii));
    % Cancel n*log(r) against n*log(beta) analytically, not numerically.
    % The (y-n)*log(beta) term MUST remain for legacy fractional inputs.
    logTerms(ii) = ni.*logMu(ii)+(y(ii)-ni).*logBeta(ii)+ ...
        correction-y(ii).*logOnePlusBeta(ii)-penalty(ii);
end
if all(isfinite(logTerms))
    nll = -sum(logTerms);
    if ~isfinite(nll), nll = Inf; end
end
end

function correction = scaledLogRisingFactorial(n,logR)
% log Gamma(r+n)-log Gamma(r)-n*log(r), for r>=16 and integer n>=0.
% Write t=n/r and H(t)=((1+t)*log1p(t)-t)/t. Then the Stirling main term is
% n*H(t)-0.5*log1p(t). This avoids subtracting quantities of order r*log(r).
% Corrections through r^(-9) have first omitted term <1.1e-16 at r=16.
correction = zeros(size(n));
active = n>1; % n=0 and n=1 give exactly zero, including at r=Inf.
if ~any(active), return; end
nn = n(active);
u = exp(-logR(active));
t = nn.*u;
L = log1p(t);
H = zeros(size(t));
small = t<=0.01;
s = t(small);
% Taylor polynomial through t^9; first omitted term <= t^10/110.
H(small) = s.*(1/2+s.*(-1/6+s.*(1/12+s.*(-1/20+s.*(1/30+ ...
    s.*(-1/42+s.*(1/56+s.*(-1/72+s./90))))))));
H(~small) = (1+1./t(~small)).*L(~small)-1;
% S(r+n)-S(r); expm1 retains accuracy when n/r is tiny.
u2 = u.*u;
stirlingDifference = u.*(expm1(-L)./12 + u2.*( ...
    -expm1(-3*L)./360 + u2.*(expm1(-5*L)./1260 + ...
    u2.*(-expm1(-7*L)./1680 + u2.*expm1(-9*L)./1188))));
correction(active) = nn.*H-0.5*L+stirlingDifference;
end
