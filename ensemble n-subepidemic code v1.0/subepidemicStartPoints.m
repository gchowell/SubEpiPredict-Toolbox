function starts = subepidemicStartPoints(map,nLHS,nJitter,seed)
%SUBEPIDEMICSTARTPOINTS Start points for MultiStart, one per row, original units.
%   starts = subepidemicStartPoints(map,nLHS,nJitter,seed)
%
%   map      from subepidemicUnitMap('make',...)
%   nLHS     Latin-hypercube points over the whole box, drawn in UNIT
%            coordinates, so a log-scale parameter is sampled log-uniformly
%            (lhsdesign when available, otherwise uniform random numbers)
%   nJitter  copies of SEED with every free parameter multiplied by
%            10^N(0, 0.1), about +/- 25% for one standard deviation; a free
%            parameter that is zero gets an additive perturbation of 1% of
%            its range instead
%   seed     a parameter vector (default map.z0); it is always the last
%            row. Fixed parameters keep the value in map.z0.
%
%   Duplicate rows are removed. The global random stream is used, so set the
%   seed with rng before fitting for reproducible start points.
%
%   See also subepidemicUnitMap, fitSubepidemicModel

if nargin<4 || isempty(seed)
    seed = map.z0;
end
seed = double(full(seed(:).'));
if numel(seed)~=map.nz || ~isreal(seed) || any(~isfinite(seed))
    error('SubEpiPredict:InvalidSeed', ...
        'The seed must be a finite real vector with %d entries.',map.nz);
end
seed = min(max(seed,map.LB),map.UB);
fixed = map.UB<=map.LB;
seed(fixed) = map.z0(fixed);

% Latin hypercube over the unit box of the free coordinates
S1 = zeros(nLHS,map.nz);
if nLHS>0 && map.nFree>0
    if exist('lhsdesign','file')==2
        X = lhsdesign(nLHS,map.nFree,'criterion','maximin','iterations',50);
    else
        X = rand(nLHS,map.nFree);
    end
    for i = 1:nLHS
        S1(i,:) = subepidemicUnitMap('toZ',X(i,:),map);
    end
elseif nLHS>0
    S1 = repmat(map.z0,nLHS,1);
end

% jittered copies of the seed
S2 = zeros(nJitter,map.nz);
if nJitter>0
    span = map.UB-map.LB;
    jitter = 10.^(0.1*randn(nJitter,map.nz));
    S2 = repmat(seed,nJitter,1).*jitter;
    zeroCols = (seed==0);
    if any(zeroCols)
        S2(:,zeroCols) = repmat(seed(zeroCols),nJitter,1)+ ...
            0.01*randn(nJitter,sum(zeroCols)).*repmat(span(zeroCols),nJitter,1);
    end
    S2 = min(max(S2,repmat(map.LB,nJitter,1)),repmat(map.UB,nJitter,1));
    S2(:,fixed) = repmat(map.z0(fixed),nJitter,1);
end

starts = [S1;S2;seed];
starts = unique(starts,'rows','stable');
keep = all(isfinite(starts),2);
starts = starts(keep,:);
if isempty(starts)
    starts = seed;
end
end
