function [rs,ps,as,Ks,alpha,d]=unpackSubepidemicParameters(z,npatches)
%UNPACKSUBEPIDEMICPARAMETERS Decode a dimension-checked parameter vector.
%   Layout: [r(1:n), p(1:n), a(1:n), K(1:n), alpha, d], exactly 4*n+2
%   entries, including fixed parameters. Accepts row or column vectors and
%   returns each growth block as a row. This helper does not use globals.
%   This is a dimension/type check, not a general optimizer-validity check.

if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SubEpiPredict:InvalidPatchCount', ...
        'npatches must be a positive integer when decoding parameters.');
end
expectedDimension=4*npatches+2;
if ~(isnumeric(z) && isreal(z) && isvector(z) && numel(z)==expectedDimension)
    error('SubEpiPredict:ParameterDimensionMismatch', ...
        'Expected a real parameter vector of length %d for %d subepidemics; received %d entries.', ...
        expectedDimension,npatches,numel(z));
end

z=z(:).';
rs=z(1:npatches);
ps=z(npatches+1:2*npatches);
as=z(2*npatches+1:3*npatches);
Ks=z(3*npatches+1:4*npatches);
alpha=z(end-1);
d=z(end);
end
