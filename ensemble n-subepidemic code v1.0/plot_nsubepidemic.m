% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function totinc=plot_nsubepidemic(flag1_pass,r_pass,p_pass,a_pass,K_pass,npatches_pass,onset_thr_pass,I0_pass,windowsize1_pass)

% Plot model simulation with a specific set of parameter values provided by
% the user

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 %Parameter estimation method - LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3,MLE (Neg Binomial)=4, MLE (Neg Binomial)=5
global npatches_fixed
global onset_fixed
global dist1
global smoothfactor1
global calibrationperiod1
global invasions

global invasions
global timeinvasions
global Cinvasions

% <============================================================================>
% <=================== Parameter values =======================================>
% <============================================================================>


if exist('flag1_pass','var')==1 & isempty(flag1_pass)==0
    flag1=flag1_pass;
else
    'flag1 value is required'
    return
end

if exist('npatches_pass','var')==1 & isempty(npatches_pass)==0
    npatches=npatches_pass;
else
    'npatches value is required'
end

if exist('onset_thr_pass','var')==1 & isempty(onset_thr_pass)==0
    onset_thr=onset_thr_pass;
else
    'onset_thr value is required'
end

% <==================================================================================>
% <========================== Parameter values ======================================>
% <==================================================================================>


a1=ones(npatches,1);
p1=ones(npatches,1);
K1=ones(npatches,1);

switch flag1

    case 1

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required'
            return
        end

        if exist('p_pass','var')==1 & isempty(p_pass)==0
            p1=p_pass;
        else
            'p value is required'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required'
            return

        end

        model_name1='Generalized Logistic Growth Model';

        params=[r1 p1 K1];

    case 2

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required'
            return
        end

        if exist('p_pass','var')==1 & isempty(p_pass)==0
            p1=p_pass;
        else
            'p value is required'
            return
        end

        if exist('a_pass','var')==1 & isempty(a_pass)==0
            a1=a_pass;
        else
            'a value is required'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required'
            return
        end

        if exist('q_pass','var')==1 & isempty(q_pass)==0
            q1=q_pass;
        else
            'q value is required'
            return
        end

        model_name1='Generalized Richards Model';

        params=[r1 p1 a1 K1];

    case 3

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required'
            return

        end

        if exist('q_pass','var')==1 & isempty(q_pass)==0
            q1=q_pass;
        else
            'q value is required'
            return
        end

        model_name1='Logistic Growth Model';

        params=[r1 K1];

    case 4

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required'
            return
        end

        if exist('a_pass','var')==1 & isempty(a_pass)==0
            a1=a_pass;
        else
            'a value is required'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required'
            return

        end

        if exist('q_pass','var')==1 & isempty(q_pass)==0
            q1=q_pass;
        else
            'q value is required'
            return
        end

        model_name1='Richards Model';

        params=[r1 a1 K1];

end

color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];

% <==================================================================================>
% <========================== Initial condition C(0) ================================>
% <==================================================================================>

if exist('I0_pass','var')==1 & isempty(I0_pass)==0
    I0=I0_pass;
else
    'Initial value C(0) is required'
    return
end

% <==================================================================================>
% <========================== Time span =============================================>
% <==================================================================================>

if exist('windowsize1_pass','var')==1 & isempty('windowsize1_pass')==0
    windowsize1=windowsize1_pass;
else
    'Simulation duration is required'
end

timevect=0:1:windowsize1;

IC=zeros(npatches,1);

if onset_thr>0
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


[~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],r1,p1,a1,K1,npatches,onset_thr,flag1);

for j=1:npatches

    incidence1=[x(1,j);diff(x(:,j))];

    line1=plot(timevect,incidence1,color1(j,:))

    set(line1,'Linewidth',2)
    hold on

end

y=sum(x,2);
totinc=[y(1,1);diff(y(:,1))];
if onset_thr>0
    totinc(1)=totinc(1)-(npatches-1);
end
bestfit=totinc;

line1=plot(timevect,totinc,'k')
set(line1,'Linewidth',3)

ylabel('dC(t)/dt')
xlabel('Time');

set(gca,'FontSize',GetAdjustedFontSize)
set(gcf,'color','white')

axis([0 timevect(end) 0 max(totinc)+10])
