function [P,group,changed] = subepidemicOrderLabels(P,npatches,onset_fixed,onset_thr,I0)
%SUBEPIDEMICORDERLABELS Fixed order for sub-epidemics that start together.
%   [P,group,changed] = subepidemicOrderLabels(P,npatches,onset_fixed,onset_thr,I0)
%
%   P        parameter vectors, one per row: [r(1:n) p(1:n) a(1:n) K(1:n) alpha d]
%   group    the sub-epidemics whose labels are interchangeable (empty: none)
%   changed  for each row, true if its labels were reordered
%
%   WHY. Sub-epidemics that start at time zero from the same value are
%   interchangeable: swapping their parameters gives exactly the same curve
%   and the same likelihood. The optimizer then returns them in either
%   order. The ranked fit and each bootstrap refit can differ in that order,
%   which mixes the two sub-epidemics in every interval of r, p, a and K and
%   makes "K_1" of one fit the "K_2" of another.
%
%   WHEN labels are interchangeable (simulateSubepidemic.m,
%   initializeSubepidemicState.m):
%     onset_fixed = 1 (synchronous)   all sub-epidemics: each starts at I0/n
%     onset_fixed = 0 and the onset threshold is already met at time zero
%       (onset_thr <= I0 for 2 sub-epidemics, <= min(I0,1) for more):
%       the first starts at I0 and every later one at 1, so
%         I0 = 1       all sub-epidemics
%         otherwise    sub-epidemics 2 to n (3 or more sub-epidemics only)
%   In every other case the order is set by the onsets and nothing is done.
%
%   RULE. Within the group the sub-epidemics are put in increasing order of
%   r (ties: increasing K). The curve, the likelihood and the AICc do not
%   change. Rows with a missing value in the group are left as they are.
%
%   Added in revision m (2026-11-18).

npatches = full(double(npatches));
nz = 4*npatches+2;
wasSparse = issparse(P);
P = full(double(P));
if isvector(P)
    P = P(:).';
end
if size(P,2)~=nz
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Each row must hold %d parameters; received %d columns.',nz,size(P,2));
end
changed = false(size(P,1),1);
group = [];
if npatches>=2
    if onset_fixed==1
        group = 1:npatches;
    else
        if npatches==2
            levelAtStart = I0;
        else
            levelAtStart = min(I0,1);
        end
        if onset_thr<=levelAtStart
            if I0==1
                group = 1:npatches;
            elseif npatches>=3
                group = 2:npatches;
            end
        end
    end
end
if numel(group)>=2
    rCols = group;
    KCols = 3*npatches+group;
    for i = 1:size(P,1)
        if any(~isfinite(P(i,[rCols KCols])))
            continue
        end
        [~,order] = sortrows([P(i,rCols)' P(i,KCols)']);
        order = order';
        if ~isequal(order,1:numel(group))
            for block = 0:3
                P(i,block*npatches+group) = P(i,block*npatches+group(order));
            end
            changed(i) = true;
        end
    end
end
if wasSparse
    P = sparse(P);
end
end
