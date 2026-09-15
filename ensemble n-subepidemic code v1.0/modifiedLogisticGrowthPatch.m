function dx=modifiedLogisticGrowthPatch(t,x,rs1,ps1,as1,Ks1,npatches,onset_thr,flag1,activeMask)
%MODIFIEDLOGISTICGROWTHPATCH Pure derivative for one fixed-activation segment.
%   ACTIVEMASK is an explicit, immutable input during each ode15s solve.
%   This function neither reads nor changes any global activation state.
%   Threshold decisions and event restarts belong to simulateSubepidemic.
%   T and ONSET_THR remain in the interface, but do not determine activation.
%   The six growth expressions are unchanged from the supplied toolbox.
%
%   Direct legacy nine-argument calls are intentionally rejected. Use
%   simulateSubepidemic for a complete asynchronous trajectory, or supply an
%   explicit tenth argument for a segment with a fixed activation mask.

if nargin<10
    error('SubEpiPredict:MissingActivationMask', ...
        'Supply an explicit activeMask, or call simulateSubepidemic.');
end
if ~((isnumeric(activeMask) || islogical(activeMask)) && ...
        isreal(activeMask) && isvector(activeMask) && ...
        numel(activeMask)==npatches && all(activeMask(:)==0 | activeMask(:)==1))
    error('SubEpiPredict:InvalidActivationMask', ...
        'activeMask must contain one zero or one for each patch.');
end

dx=zeros(npatches,1);
for j=1:npatches
    if ~activeMask(j)
        % Do not evaluate a dormant growth expression: 0*NaN is still NaN.
        continue
    end
    switch flag1
        case 0 % GGM
            dx(j)=rs1(j)*x(j).^ps1(j);
        case 1 % GLM
            dx(j)=(rs1(j)*x(j).^ps1(j))*(1-x(j)/Ks1(j));
        case 2 % Generalized Richards (existing expression retained)
            dx(j)=(rs1(j)*x(j).^ps1(j))*(1-x(j)/Ks1(j)).^as1(j);
        case 3 % Logistic
            dx(j)=(rs1(j)*x(j))*(1-x(j)/Ks1(j));
        case 4 % Richards
            dx(j)=(rs1(j)*x(j))*(1-(x(j)/Ks1(j)).^as1(j));
        case 5 % Gompertz
            dx(j)=rs1(j)*x(j)*log(Ks1(j)/x(j));
        otherwise
            error('SubEpiPredict:InvalidModelFlag','flag1 must be an integer from 0 to 5.');
    end
end
end
