function [AICc,part1,part2,numparams]=getAICc(method1,dist1,npatches,flag1,fixI0,fval,n,onset_fixed)


% last revised: 01 June 2022

%flag 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz

numparams=0;

for j=1:npatches

    switch flag1 % model indicator

        case 0
            numparams=numparams+2;

        case 1
            numparams=numparams+3;

        case 2
            numparams=numparams+4;

        case 3
            numparams=numparams+2;

        case 4
            numparams=numparams+3;

        case 5
            numparams=numparams+2;

    end
end


if npatches>1 & onset_fixed==0

    numparams=numparams+1; %Cthr

end

if fixI0==0 %fix initial datum or estimated

    numparams=numparams+1;

end


if method1==0 & dist1==0 % Normal distribution -- one parameter for variance

    numparams=numparams+1;

elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)

    numparams=numparams+1;

elseif method1==5

    numparams=numparams+2;  %Neg. Binomial requires 2 more parameters (alpha,d)

elseif method1==6  % Laplace distribution -- one parameter for scale

    numparams=numparams+1;

end


% Revision h. The small-sample term 2k(k+1)/(n-k-1) is defined only for
% n > k+1. At n = k+1 it is infinite and below that it is NEGATIVE, so a model
% with more parameters than the data can support was rewarded (with 10 points
% and least squares: +16 for one GLM sub-epidemic, -110 for three). Such a
% model now gets AICc = Inf and can never be ranked ahead of a valid one. The
% screen (fittingModifiedLogisticFunctionPatchABC.m) does not fit it at all.
if n-numparams-1<=0
    AICc=Inf;
    part1=NaN;
    part2=Inf;
    return
end

switch method1

    case 0

        AICc= n*log(fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

        part1= n*log(fval);

        part2= 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

    case 6

        AICc= 2*n*log(2*fval/n) + 2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

        part1 = 2*n*log(2*fval/n);

        part2=  2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);


    otherwise

        AICc=-2*(-fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1=-2*(-fval);
        
        part2=2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

end


%method1

%fval

%n

%numparams
