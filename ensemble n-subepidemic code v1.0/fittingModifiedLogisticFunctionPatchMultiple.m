function [Phatss,npatches,onset_thr,curves,bestfit,data1,P0,AICc_best,factor1,d]=fittingModifiedLogisticFunctionMultiple(RMSES,PS,data1,DT,epidemic_period,M,flagX,numstartpoints,rank1)

global flag1 timevect ydata yfit calibrationperiod1

global I0 npatches onset_thr

flag1=flagX;

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global method1
global dist1
global factor1

global smoothfactor1

global LBe UBe


close all

% <=============================================================================================>
% <========== Load previously inferred models according to AICc (best to worst fits) ===========>
% <=============================================================================================>

%load(strcat('./output/ABC-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4),'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

% remove repeated rows
[RMSES,index1]=unique(RMSES,'rows','stable');
PS=PS(index1,:);

index1=rank1;

%[npatches,onset_thr,P0]

npatches=RMSES(index1,1);

onset_thr=RMSES(index1,2);

AICc_best=RMSES(index1,3);

P0=PS(index1,1:npatches*4+2);

numparams=get_nparams(method1,npatches,flag1,1);

%

% fitting

Phats=[];

goodnessfit=[];

modelfits=[];

SSstats1=[];

SSstats2=[];

%data1=load(strcat('./output/',datafilename1));

data1=data1(epidemic_period,:);

data1(:,2)=data1(:,2);

I0=data1(1,2); % initial condition

if I0==0
    
    data1=data1(2:end,:);
    
end

data=data1(:,2);

[max1,index1]=max(data);

timevect=(data1(:,1))*DT;

rs1=P0(1:npatches);

ps1=P0(npatches+1:2*npatches);

as1=P0(2*npatches+1:3*npatches);

Ks1=P0(3*npatches+1:4*npatches);

alpha=P0(end-1);

d=P0(end);


z=[rs1 ps1 as1 Ks1 alpha d];

[LB1,UB1]=getbounds(npatches);

LB=[LB1 LBe]; %r p a K alpha d

UB=[UB1 UBe];



I0=data(1); % initial condition

hold on

%options=optimset('MaxFunEvals',3000,'MaxIter',3000,'Algorithm','trust-region-reflective','TolFun',1e-6,'TolX',1e-6);

%options=[];

%options = optimoptions('lsqcurvefit','UseParallel',true,...
%    'TolX',10^(-5),'TolFun',10^(-5));

%[P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowthPatch1,z,timevect,smooth(data,smoothfactor1),LB,UB,options,I0,npatches,onset_thr,flag1);

% resnorm is the SSE which is given by sum(residual.^2)
% P is the vector with the estimated parameters


%method1=3; %LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3


ydata=smooth(data,smoothfactor1);

options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);

%options=optimoptions('fmincon','Algorithm','sqp','tolfun',10^-6,'TolX',10^-6,'MaxFunEvals',3200,'MaxIter',3200);

f=@plotModifiedLogisticGrowthPatchMethodsFullLogLik;

problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);

%ms = MultiStart('PlotFcns',@gsplotbestf);
ms = MultiStart('Display','final');
%ms=MultiStart;

%rpoints = RandomStartPointSet('NumStartPoints',numstartpoints); % start with a few random starting sets in addition to the guess supplied by the user (z)

pts = z;
tpoints = CustomStartPointSet(z);
%allpts = {tpoints,rpoints};
allpts = {tpoints};

[P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);


rs_hat=P(1,1:npatches);
ps_hat=P(1,npatches+1:2*npatches);
as_hat=P(1,2*npatches+1:3*npatches);
Ks_hat=P(1,3*npatches+1:4*npatches);

alpha_hat=P(end-1);
d_hat=P(end);

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


figure(10)

for j=1:npatches
    
    incidence1=[x(1,j);diff(x(:,j))];
    
    plot(timevect,incidence1)
    hold on
    
end

y=sum(x,2);

totinc=[y(1,1);diff(y(:,1))];

if onset_fixed==0
    totinc(1)=totinc(1)-(npatches-1);
end

bestfit=totinc;

% normal distribution of the error structure
if method1==0 & dist1==0

    var1=sum((bestfit-data).^2)./(length(bestfit)-numparams); % last revised: 01 June 2022
    factor1=sqrt(var1);

end

        
plot(timevect,totinc,'r')
hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('incidence fit')

title(strcat('Num. Subepidemics=',num2str(npatches),'; AICc=',num2str(AICc_best,6)))

legend(strcat('Sub-epidemics=',num2str(npatches),'; C_{thr}=',num2str(onset_thr)))
set(gca,'FontSize',24)
set(gcf,'color','white')


% <=============================================================================================>
% <======================== Vector with parameter estimates ====================================>
% <=============================================================================================>

Ptrue=[rs_hat ps_hat as_hat Ks_hat alpha_hat d_hat]

'generate simulation study to derive parameter uncertainty..'

% <=============================================================================================>
% <========== Generate parameter uncertainty via parametric bootstrapping ======================>
% <=============================================================================================>

yi=cumsum(totinc);

z=Ptrue; %r p a Ks alpha(MLE-Neg binomial) d

Phatss=zeros(M,4*npatches+2);
SSEs=zeros(M,1);
curves=zeros(length(yi),M);


% <================================================================================================>
% <=========================== Set initial parameter guesses and bounds ===========================>
% <================================================================================================>

rs1=P0(1:npatches);

ps1=P0(npatches+1:2*npatches);

as1=P0(2*npatches+1:3*npatches);

Ks1=P0(3*npatches+1:4*npatches);

alpha=P0(end-1);

d=P0(end);

z=[rs1 ps1 as1 Ks1 alpha d];

[LB1,UB1]=getbounds(npatches);

LB=[LB1 LBe]; %r p a K alpha d

UB=[UB1 UBe];
 

for real=1:M
     
    yirData=AddPoissonError(yi,1,dist1,factor1,d_hat);
    
    curves(:,real)=yirData;
    
    
    I0=data(1); % initial condition
    
    ydata=yirData;

    %ydata=smooth(yirData,smoothfactor1);
    
    %options=optimset('MaxFunEvals',3000,'MaxIter',3000,'Algorithm','trust-region-reflective','TolFun',1e-5,'TolX',1e-5);
    
    %options=[];
    
    %options = optimoptions('lsqcurvefit','UseParallel',true,...
    %    'TolX',10^(-5),'TolFun',10^(-5));
    
    options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);
    
    f=@plotModifiedLogisticGrowthPatchMethodsFullLogLik;
    
    problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);
    
    %ms = MultiStart('PlotFcns',@gsplotbestf);
    %ms = MultiStart('Display','final');
    ms = MultiStart('Display','off');
    
    %ms=MultiStart;
    
    %one random start guess
    rpoints = RandomStartPointSet('NumStartPoints',2); % start with a few random starting sets in addition to the guess supplied by the user (z)
    
    pts = z;
    tpoints = CustomStartPointSet(z);
    allpts = {tpoints,rpoints};
    %allpts = {tpoints};
    
    [P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);
    
    Phatss(real,:)=P;
    
end

