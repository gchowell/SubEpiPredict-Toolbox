% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [RMSES,PS,npatches,onset_thr,P]=fittingModifiedLogisticFunctionPatchABC(datafilename1,data1,DT,epidemic_period,M,flagX,numstartpoints)

global flag1 method1 timevect ydata

global I0 npatches onset_thr

flag1=flagX;

close all

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global dist1
global factor1

global smoothfactor1
global calibrationperiod1

global LBe UBe

% <=============================================================================================>
% <========= Set bounds for parameters associated with error structure (alpha and d) ===========>
% <=============================================================================================>

switch method1
    
    case 0
        LBe=[0 0];
        UBe=[0 0];
    case 1
        LBe=[0 0];
        UBe=[0 0];
    case 2
        LBe=[0 0];
        UBe=[0 0];
    case 3
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 4
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 5
        LBe=[10^-8 0.2];
        UBe=[10^5 10^2];
end

% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>

%data1=load(strcat('./output/',datafilename1));

data1=data1(epidemic_period,:);

data1(:,2)=data1(:,2);

I0=data1(1,2); % initial condition

if I0==0
    data1=data1(2:end,:);
end

data=data1(:,2);

% <==============================================================================>
% <============ Set time vector (timevect) and initial condition (I0) ===========>
% <==============================================================================>

timevect=(data1(:,1));

I0=data(1); % initial condition

% <==============================================================================>
% <===================== Set initial parameter guesses ==========================>
% <==============================================================================>

rs1=zeros(1,npatches_fixed)+0.1;
ps1=zeros(1,npatches_fixed)+0.9;
Ks1=zeros(1,npatches_fixed)+sum(data1(:,2));

as1=ones(1,npatches_fixed);

for j=1:npatches_fixed
    
    if flag1==3 | flag1==4 | flag1==5 % Logistic model or Richards model (p=1)
        ps1(j)=1;
    else
        ps1(j)=0.9;
    end
    
    if flag1==5
        
        rs1(j)=1-I0/Ks1(j);
        
        as1(j)=rs1(j)/log(Ks1(j)/I0);
        
    end
    
end

% <==============================================================================>
% <================= Set range of C_thr values (onset_thrs) =====================>
% <==============================================================================>

%cum1=sum(smooth(data1(:,2),smoothfactor1));
%onset_thrs=unique(cumsum(smooth(data1(:,2),smoothfactor1)));
%index2=find(onset_thrs<=0.99*cum1);
%onset_thrs=onset_thrs(index2)';

% equal-witdh discretization of C_thr
cumcurve1=cumsum(smooth(data1(:,2),smoothfactor1));

onset_thrs=linspace(cumcurve1(1),cumcurve1(end),length(data1(:,2)));

onset_thrs=[0 onset_thrs(1:end-1)];

% <=============================================================================>
% <===== Set range of the possible number of subepidemics (1:npatches_fixed)====>
% <=============================================================================>

npatchess=1:1:npatches_fixed;

if onset_fixed==1 | (length(npatchess)==1 & npatchess(1)==1)
    
    onset_thrs=0;
    
end

onset_thrs2=onset_thrs;

RMSES=sparse(1000,3);

PS=sparse(1000,npatches_fixed*4+2);

count1=1;

% <====================================================================================>
% <==== Evaluate AICc across models with different number of subepidemics and C_thr ==>
% <====================================================================================>

ydata=smooth(data,smoothfactor1);

for npatches2=[npatchess]
    
    npatches=npatches2;
    
    if (onset_fixed==1 | npatches==1)
        onset_thrs=0;
    else
        onset_thrs=onset_thrs2;
    end
    
    % <================================================================================================>
    % <=========================== Set initial parameter guesses and bounds ===========================>
    % <================================================================================================>
    
    rs1=zeros(1,npatches)+0.1;
    ps1=zeros(1,npatches)+0.9;
    Ks1=zeros(1,npatches)+sum(data1(:,2));
    
    as1=ones(1,npatches);
    
    for j=1:npatches
        if flag1==3 | flag1==4 | flag1==5 % Logistic model or Richards model (p=1)
            ps1(j)=1;
        else
            ps1(j)=0.9;
        end
        if flag1==5
            rs1(j)=1-I0/Ks1(j);
            as1(j)=rs1(j)/log(Ks1(j)/I0);
        end
    end
    
    z=[rs1 ps1 as1 Ks1 1 1];
    
    [LB1,UB1]=getbounds(npatches,data);
    
    LB=[LB1 LBe]; %r p a K alpha d
    UB=[UB1 UBe];
    
    % --- safety: flip any LB>UB just in case --------------------------------
    flipMask = LB > UB;
    if any(flipMask)
        tmp = LB(flipMask); LB(flipMask) = UB(flipMask); UB(flipMask) = tmp;
    end
    % ------------------------------------------------------------------------
    
    % ========================= CHANGED [2] START =============================
    % Generate MultiStart initial points ONCE per (npatches, LB/UB)
    z0   = min(max(z,LB),UB);                 % clamp seed into bounds        %CHANGED [2]
    dDim = numel(LB);                                                      %CHANGED [2]
    span = (UB - LB);                                                     %CHANGED [2]
    nStarts = max(0, floor(numstartpoints));                              %CHANGED [2]
    
    if nStarts > 0                                                         %CHANGED [2]
        if exist('lhsdesign','file') == 2                                  %CHANGED [2]
            X = lhsdesign(nStarts, dDim, 'criterion','maximin','iterations',50); %CHANGED [2]
        else                                                               %CHANGED [2]
            X = rand(nStarts, dDim);                                       %CHANGED [2]
        end                                                                %CHANGED [2]
        starts_base = LB + X .* span;                                      %CHANGED [2]
        starts_base = [starts_base; z0];    % include user seed             %CHANGED [2]
    else                                                                   %CHANGED [2]
        starts_base = z0;                                                  %CHANGED [2]
    end                                                                    %CHANGED [2]
    starts_base = unique(round(starts_base,6), 'rows');                    %CHANGED [2]
    inB = all(starts_base >= (LB - 1e-12) & starts_base <= (UB + 1e-12), 2); %CHANGED [2]
    finiteReal = all(isfinite(starts_base),2) & isreal(starts_base);       %CHANGED [2]
    starts_base = starts_base(inB & finiteReal, :);                        %CHANGED [2]
    if isempty(starts_base), starts_base = z0; end                         %CHANGED [2]
    % ========================== CHANGED [2] END ==============================
    
    nloops=length(onset_thrs);
    
    RMSES2=sparse(1000,3);
    count2=1;

    for onset_thr=onset_thrs
        % ******** MLE estimation method with MultiStart  *********
        % check multiple initial guesses to ensure global minimum is obtained
               
        options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);
        
        f=@plotModifiedLogisticGrowthPatchMethodLogLik;
        
        % ========================= CHANGED [2] START =========================
        % Reuse the base start set for every onset_thr                        %CHANGED [2]
        problem = createOptimProblem('fmincon','objective',f,'x0',z0,'lb',LB,'ub',UB,'options',options); %CHANGED [2]
        useParallel = ~isempty(gcp('nocreate'));                             %CHANGED [2]
        ms = MultiStart('Display','off','UseParallel',useParallel,'StartPointsToRun','bounds-ineqs'); %CHANGED [2]
        sp = CustomStartPointSet(starts_base);                               %CHANGED [2]
        % ========================== CHANGED [2] END ==========================
        
        [P,fval,flagg,outpt,allmins] = run(ms,problem,sp);
        
        % --> numerical solver to get the best fit in order to check the actual number of
        % subepidemics involved in the best fit
        
        rs_hat=P(1,1:npatches);
        ps_hat=P(1,npatches+1:2*npatches);
        as_hat=P(1,2*npatches+1:3*npatches);
        Ks_hat=P(1,3*npatches+1:4*npatches);
        
        alpha_hat=P(1,end-1);
        d_hat=P(1,end);
        
        IC=zeros(npatches,1);
        
        if onset_fixed==0
            IC(1,1)=I0;
            IC(2:end,1)=1;
            
            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);
            
            invasions(1)=1;
            timeinvasions(1)=0;
            Cinvasions(1)=0;
        else
            IC(1:end,1)=I0./length(IC(1:end,1));
            
            invasions=zeros(npatches,1);
            timeinvasions=zeros(npatches,1);
            Cinvasions=zeros(npatches,1);
            
            invasions(1:end)=1;
            timeinvasions(1:end)=0;
            Cinvasions(1:end)=0;
        end
        
        [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);
        
        if sum(invasions)==1 & sum(invasions)<npatches
            continue
        elseif sum(invasions)<npatches
            npatches=sum(invasions);
            P=[rs_hat(1:npatches) ps_hat(1:npatches) as_hat(1:npatches) Ks_hat(1:npatches) alpha_hat d_hat];
            %pause
        end
        
        AICc=getAICc(method1,dist1,npatches,flag1,1,fval,length(ydata),onset_fixed);
        
        RMSES(count1,:)=[npatches onset_thr AICc];
        PS(count1,1:length(P))=P;
        count1=count1+1;
                
    end %onset
    
end %npatches

%RMSES(1:count1,:)
%pause

% <=============================================================================================>
% <======================== Sort the results by AICc (lowest to highest) =======================>
% <=============================================================================================>

RMSES=RMSES(1:count1-1,:);

PS=PS(1:count1-1,:);

[RMSES,index1]=sortrows(RMSES,[3 1]);

PS=PS(index1,:);

[RMSE1, index1]=min(RMSES(:,3));

npatches=RMSES(index1,1);

onset_thr=RMSES(index1,2);

AICc_best=RMSES(index1,3);

%-->If we have a series of AICc  (or wSSE) values from N different models sorted from lowest (best model) 
%to highest (worst model), I am wondering if we could define a proper threshold criterion to drop models with associated AICc (or wSSE) value greater than some threshold criteria.

% -->Let AICmin denote the minimun AIC from several models. 
%
