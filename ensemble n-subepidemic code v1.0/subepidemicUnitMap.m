function out = subepidemicUnitMap(action,arg,map)
%SUBEPIDEMICUNITMAP Map between parameter vectors and unit coordinates.
%   The optimizer does not see the full vector z = [r, p, a, K, alpha, d].
%   It sees only the FREE entries (UB > LB), each mapped to [0, 1]:
%       linear coordinate   u = (z - LB) / (UB - LB)
%       log coordinate      u = (log z - log LB) / (log UB - log LB)
%   so every coordinate has the same scale whatever the size of the
%   parameter. One function holds both directions, so the fit, the start
%   points and the diagnostics cannot disagree about the mapping.
%
%   map = subepidemicUnitMap('make',struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog))
%       builds the map. z0 is projected into the bounds; its fixed entries
%       are the values every mapped vector keeps.
%   U = subepidemicUnitMap('toU',Z,map)
%       Z holds one parameter vector per ROW (original units); U holds the
%       unit coordinates of the free entries, one point per row, clipped to
%       [0, 1].
%   z = subepidemicUnitMap('toZ',u,map)
%       u is ONE point in unit coordinates; z is the full row vector in
%       original units, inside the bounds.
%
%   See also fitSubepidemicModel, subepidemicStartPoints, subepidemicSearchSpace

switch action
    case 'make'
        z0 = double(full(arg.z0(:).'));
        LB = double(full(arg.LB(:).'));
        UB = double(full(arg.UB(:).'));
        isLog = logical(full(arg.isLog(:).'));
        nz = numel(LB);
        if any([numel(z0),numel(UB),numel(isLog)]~=nz)
            error('SubEpiPredict:ParameterDimensionMismatch', ...
                'Seed, bounds and scale flags must have the same length.');
        end
        if any(~isfinite(LB)) || any(~isfinite(UB)) || any(LB>UB)
            error('SubEpiPredict:ReversedBounds', ...
                'Bounds must be finite with LB <= UB. Reversed bounds are not swapped.');
        end
        z0 = min(max(z0,LB),UB);
        if ~isreal(z0) || any(~isfinite(z0))
            error('SubEpiPredict:InvalidSeed','The seed is not a finite real vector.');
        end
        isLog = isLog & LB>0 & UB>LB;
        gLB = LB; gUB = UB;
        gLB(isLog) = log(LB(isLog));
        gUB(isLog) = log(UB(isLog));
        out = struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog, ...
            'free',find(UB>LB),'gLB',gLB,'span',gUB-gLB,'nz',nz);
        out.nFree = numel(out.free);

    case 'toU'
        Z = double(full(arg));
        if isvector(Z), Z = Z(:).'; end
        if size(Z,2)~=map.nz
            error('SubEpiPredict:ParameterDimensionMismatch', ...
                'Each row must hold %d parameters; received %d columns.',map.nz,size(Z,2));
        end
        G = Z;
        if any(map.isLog)
            G(:,map.isLog) = log(max(Z(:,map.isLog),realmin));
        end
        n = size(Z,1);
        out = (G(:,map.free)-repmat(map.gLB(map.free),n,1))./repmat(map.span(map.free),n,1);
        out = min(max(out,0),1);

    case 'toZ'
        u = double(arg(:).');
        if numel(u)~=map.nFree
            error('SubEpiPredict:ParameterDimensionMismatch', ...
                'Expected %d unit coordinates; received %d.',map.nFree,numel(u));
        end
        u = min(max(u,0),1);
        g = map.gLB(map.free)+u.*map.span(map.free);
        lg = map.isLog(map.free);
        g(lg) = exp(g(lg));
        out = map.z0;
        out(map.free) = min(max(g,map.LB(map.free)),map.UB(map.free));

    otherwise
        error('SubEpiPredict:UnitMapAction','Action must be ''make'', ''toU'' or ''toZ''.');
end
end
