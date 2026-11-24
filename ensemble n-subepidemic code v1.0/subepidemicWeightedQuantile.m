function q = subepidemicWeightedQuantile(X,p,w)
%SUBEPIDEMICWEIGHTEDQUANTILE Quantiles of weighted curves at each time point.
%   q = subepidemicWeightedQuantile(X,p,w)
%
%   X   T-by-N, one row per time point and one column per curve
%   p   probabilities, each from 0 to 1
%   w   N nonnegative weights, one per curve; they need not sum to 1.
%       [] or omitted = every curve weighs the same.
%   q   T-by-numel(p); q(i,k) is the p(k) quantile of row i of X
%
%   Definition. At one time point, sort the values x(1) <= ... <= x(n) and
%   carry their weights along, scaled to sum to 1. The value x(j) is the
%   quantile of probability
%       P(j) = w(1) + ... + w(j-1) + w(j)/2,
%   quantiles between two values are linearly interpolated, and
%   probabilities below P(1) or above P(n) give x(1) or x(n). With equal
%   weights P(j) = (j - 0.5)/n, which is the definition of quantile() of
%   the Statistics and Machine Learning Toolbox: equal weights return the
%   values of quantile(X,p,2).
%
%   Equal values share their weights equally, so the result does not depend
%   on the order of the curves (forecast counts repeat often). Curves with
%   weight 0 and NaN values are left out. A row with no value left gives NaN.
%
%   Used for the ensembles (getensemblesubepidemicmodels.m): member k with
%   weight W(k) and M(k) curves gives all its curves, each with weight
%   W(k)/M(k). Added in the revision of 2026-11-24.
%
%   See also getensemblesubepidemicmodels, computeWIS, computeQuantiles,
%   computeforecastperformance

if nargin<3
    w=[];
end
if ~(isnumeric(X) && isreal(X) && ndims(X)==2)
    error('SubEpiPredict:WeightedQuantile','X must be a real matrix, one row per time point and one column per curve.');
end
if ~(isnumeric(p) && isreal(p) && ~isempty(p) && all(p(:)>=0 & p(:)<=1))
    error('SubEpiPredict:WeightedQuantile','The probabilities must be numbers from 0 to 1.');
end
p=double(p(:)');
[T,N]=size(X);
if isempty(w)
    w=ones(1,N);
end
if ~(isnumeric(w) && isreal(w) && numel(w)==N && all(isfinite(w(:))) && all(w(:)>=0))
    error('SubEpiPredict:WeightedQuantile','Give one finite, nonnegative weight per curve (%d curves, %d weights).',N,numel(w));
end
w=full(double(w(:)'));

q=nan(T,numel(p));
used=w>0;

for i=1:T

    x=double(X(i,:));
    ok=used & ~isnan(x);
    if ~any(ok)
        continue
    end

    [xs,order]=sort(x(ok));
    ws=w(ok);
    ws=ws(order);
    ws=ws/sum(ws);
    n=numel(xs);

    % equal values share their weights equally
    group=cumsum([1 diff(xs)~=0]);
    if group(n)<n
        groupWeight=accumarray(group(:),ws(:));
        groupSize=accumarray(group(:),1);
        ws=reshape(groupWeight(group)./groupSize(group),1,n);
    end

    % probability of each sorted value: the weight below it plus half its own
    P=cumsum(ws)-ws/2;

    for k=1:numel(p)
        j=sum(P<=p(k)); % the last value at or below p(k)
        if j==0
            q(i,k)=xs(1);
        elseif j==n
            q(i,k)=xs(n);
        else
            % P(j) <= p(k) < P(j+1), so the denominator is positive
            q(i,k)=xs(j)+(p(k)-P(j))/(P(j+1)-P(j))*(xs(j+1)-xs(j));
        end
    end

end

end
