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

timevect=(data1(:,1))*DT;

I0=data(1); % initial condition


% <==============================================================================>
% <===================== Set initial parameter guesses ==========================>
% <==============================================================================>

rs1=zeros(1,npatches_fixed)+0.1;
ps1=zeros(1,npatches_fixed)+0.9;
Ks1=zeros(1,npatches_fixed)+sum(data1(:,2));

as1=ones(1,npatches_fixed);

for j=1:npatches_fixed
    
    if flag1(j)==3 | flag1(j)==4 | flag1(j)==5 % Logistic model or Richards model (p=1)
        ps1(j)=1;
    else
        ps1(j)=0.9;
    end
    
    if flag1(j)==5
        
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

onset_thrs=onset_thrs(1:end-1);

% <==============================================================================>
% <===== Set range of the possible number of subepidemics (1:npatches_fixed)=====>
% <==============================================================================>

npatchess=1:1:npatches_fixed;

if onset_fixed==1 | (length(npatchess)==1 & npatchess(1)==1)
    
    onset_thrs=0;
    
end


onset_thrs2=onset_thrs;

RMSES=sparse(1000,3);

PS=sparse(1000,npatches_fixed*4+2);

count1=1;


% <================================================================================================>
% <==== Evaluate AICc across models with different number of subepidemics and C_thr values ========>
% <================================================================================================>

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
        
        if flag1(j)==3 | flag1(j)==4 | flag1(j)==5 % Logistic model or Richards model (p=1)
            ps1(j)=1;
            
        else
            ps1(j)=0.9;
            
        end
        
        if flag1(j)==5
            
            
            rs1(j)=1-I0/Ks1(j);
            
            as1(j)=rs1(j)/log(Ks1(j)/I0);
            
        end
        
    end
    
    z=[rs1 ps1 as1 Ks1 1 1];
    
    [LB1,UB1]=getbounds(npatches);
    
    LB=[LB1 LBe]; %r p a K alpha d
    
    UB=[UB1 UBe];
    %
    
    nloops=length(onset_thrs);
    
    RMSES2=sparse(1000,3);
    count2=1;


    for onset_thr=onset_thrs
    %parfor i=1:nloops

        %onset_thr=onset_thrs(i);
        
        %[P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowthPatch1,z,timevect,smooth(data,smoothfactor1),LB,UB,options,I0,npatches,onset_thr,flag1);
        
        % ******** MLE estimation method with MultiStart  *********
        % check multiple initial guesses to ensure global minimum is obtained
        
        %method1=3; %LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3
        
        
        %'UseParallel','always'
        options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);
        
        %options=optimoptions('fmincon','Algorithm','sqp','tolfun',10^-6,'TolX',10^-6,'MaxFunEvals',20000,'MaxIter',20000);
        
        f=@plotModifiedLogisticGrowthPatchMethodsFullLogLik;
        
        problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);
        
        %ms = MultiStart('PlotFcns',@gsplotbestf);
        
        %ms = MultiStart('Display','final');
        ms = MultiStart('Display','off');
        
        %ms=MultiStart;
        
        rpoints = RandomStartPointSet('NumStartPoints',numstartpoints); % start with a few random starting sets in addition to the guess supplied by the user (z)
        
        tpoints = CustomStartPointSet(z);
        allpts = {tpoints,rpoints};
        %allpts = {tpoints};
        
        %z
        %list(tpoints)
        
        %ms = MultiStart(ms,'StartPointsToRun','bounds')
        %[xmin,fmin,flag,outpt,allmins] = run(ms,problem,allpts);
        
        
        [P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);
        
        
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

        %
        
        
        AICc=getAICc(method1,npatches,flag1,1,fval,length(ydata));
        
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
% The quantity exp((AICmin − AICi)/2) is interpreted as the relative likelihood of model i. 
% We can set an alpha (e.g., 0.05), drop the models with exp((AICmin − AICi)/2) smaller than alpha, 
% and combine other models with weighted average, 
% where the weight is proportional to exp((AICmin − AICi)/2). This is another way to assign weights,
% compared to 1/SSE. 


AICmin=RMSES(1,3);

relativelik_i=exp((AICmin-RMSES(:,3))/2);

% -> Drop models with alpha below 0.05
%index2=find(relativelik_i>0.05);
%RMSES=RMSES(index2,:);
%relativelik_i=relativelik_i(index2);

if 0

    subplot(1,2,1)
    line1=plot(RMSES(:,3),'ko-')

    set(line1,'LineWidth',2)

    xlabel('Model i')
    ylabel('AICc')
    set(gca,'FontSize', 16);
    set(gcf,'color','white')

    subplot(1,2,2)


    line1=plot(relativelik_i,'ko-')

    set(line1,'LineWidth',2)

    xlabel('Model i')
    ylabel('Relative likelihood')
    set(gca,'FontSize', 16);
    set(gcf,'color','white')

end



% <=============================================================================================>
% <============================ Get the best fit resultls ======================================>
% <=============================================================================================>


P=PS(index1,1:npatches*4+2);

rs_hat=P(1,1:npatches);
ps_hat=P(1,npatches+1:2*npatches);
as_hat=P(1,2*npatches+1:3*npatches);
Ks_hat=P(1,3*npatches+1:4*npatches);

alpha_hat=P(1,end-1);
d_hat=P(1,end);

if method1==3
    
    dist1=3; % VAR=mean+alpha*mean;
    
    factor1=alpha_hat;
    
elseif method1==4
    
    dist1=4; % VAR=mean+alpha*mean^2;
    
    factor1=alpha_hat;
    
elseif method1==5
    
    dist1=5; % VAR=mean+alpha*mean^2;
    
    factor1=alpha_hat;
    
end


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

if sum(invasions)<npatches
    
    npatches=sum(invasions);
    
    P=[rs_hat(1:npatches) ps_hat(1:npatches) as_hat(1:npatches) Ks_hat(1:npatches) alpha_hat d_hat];
    
    PS(1,:)=0;
    
    PS(1,1:length(P))=P;
    
    RMSES(1,1)=npatches;
    
end


% <=============================================================================================>
% <================================= Plot best model fit ======================================>
% <=============================================================================================>

figure(100)

for j=1:npatches
    
    incidence1=[x(1,j);diff(x(:,j))];
    
    plot(timevect,incidence1)
    hold on
    
end

y=sum(x,2);

totinc=[y(1,1);diff(y(:,1))];

if onset_thr>0
    totinc(1)=totinc(1)-(npatches-1);
end

bestfit=totinc;

plot(timevect,totinc,'r')

hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('Data and best SIR fit')

title('best fit')
[npatches onset_thr]


if method1==0
    
    if dist1==2  % calculate the overdispersion factor
        
        %     [coef,ns]=getMeanVarLinear(data,totinc,6);
        %
        %     if coef>0
        %         factor1=coef;
        %     else
        
        
        % estimate dispersion in data
        binsize1=7; %4
        
        [ratios,~]=getMeanVarianceRatio(data,binsize1,2);  % **
        
        %[ratios,~]=getMeanVarianceRatio(data,binsize1,1);
        
        index1=find(ratios(:,1)>0);
        
        factor1=mean(ratios(index1,1));
        
        factor1
        
    end
    
end


'ABC estimates:'
'[npatches onset_thr q]'

[npatches onset_thr]

% <=============================================================================================>
% <===================================  Save the results  ======================================>
% <=============================================================================================>

save(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4),'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

