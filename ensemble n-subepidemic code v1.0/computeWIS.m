function [WISC,WISFS,quantilesc1,quantilesf1]=computeWIS(data1,datalatest,curvesforecasts2,forecastingperiod,weights2)

% weights2 (optional, revision of 2026-11-24): one weight per column of
% curvesforecasts2, for the ensembles (subepidemicWeightedQuantile.m).
% Omitted or []: every curve weighs the same and quantile() is used, as before.
if nargin<5
    weights2=[];
end

%What is Alpha quantile?
%More precisely, an α-quantile is a real number Xα such that the random variable is less or equal to it with probability at least α, and greater or equal to it with probability at least (1–α). It is not possible to require equality because the probability of the value Xα may be positive.May 18, 2011

%weighted interval score

% calibration period

tf2=length(data1(:,1));

dataCperiod=1:tf2;

datac=datalatest(dataCperiod,2);

alphas=[0.02 0.05 0.1:0.1:0.9]; % consistent with COVID forecasting hub paper

alphaquantiles=[0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990];

w0=1/2;

K=length(alphas);

WISC=zeros(length(datac),1);

quantilesc1=zeros(length(datac),length(alphaquantiles));

for i=1:length(datac)
    
    sum1=0;
    
    y=datac(i);
    
    IS=zeros(length(alphas),1);
    
    for k=1:K
        
        alpha=alphas(k);
        
        w_k=alpha/2;
        
        [alpha/2  1-alpha/2];
        
        
        % (1-lpha)x100 % prediction coverage
        
        Lt=localRowQuantile(curvesforecasts2(i,:),alpha/2,weights2);

        Lt=(Lt>=0).*Lt;
        
        Ut=localRowQuantile(curvesforecasts2(i,:),1-alpha/2,weights2);

        Ut=(Ut>=0).*Ut;

        coverage1=find(y>=Lt & y<=Ut);
        
        IS(k)=[(Ut-Lt)+(2/alpha).*(Lt-y).*(y<Lt)+(2/alpha)*(y-Ut).*(y>Ut)];
        
        sum1=sum1+w_k*IS(k);
        
        
    end
    
    %median prediction, m
    m=localRowQuantile(curvesforecasts2(i,:),0.5,weights2);
    
    WISC(i)=(1/(K+1/2))*(w0*abs(y-m) + sum1);
    
    if isempty(weights2)
        quantilesc1(i,:)=quantile(curvesforecasts2(i,:),alphaquantiles);
    else
        quantilesc1(i,:)=subepidemicWeightedQuantile(curvesforecasts2(i,:),alphaquantiles,weights2);
    end
    quantilesc1(i,:)=(quantilesc1(i,:)>=0).*quantilesc1(i,:);

end %y

WISC=mean(WISC);


%Forecasting period

WISFS=[];

quantilesf1=zeros(forecastingperiod,length(alphaquantiles));

if forecastingperiod>0
    
    for forecastingperiod1=1:1:forecastingperiod
        
        dataFperiod=tf2+1:tf2+forecastingperiod1;
       
        
        datagroundtruth=datalatest;
        
        dataf=datagroundtruth(dataFperiod,2);
        
        %         if forecastingperiod1==forecastingperiod
        %             close all
        %
        %             figure
        %             plot(1:length(dataFperiod),curvesforecasts2(dataFperiod,:)','c-')
        %             hold on
        %             plot(dataf,'ro')
        %
        %             pause
        %         end
        
        alphas=[0.02 0.05 0.1:0.1:0.9];
        
        w0=1/2;
        
        WISF=zeros(length(dataf),1);
        
        for i=1:length(dataf)
            
            sum1=0;
            
            y=dataf(i);
            
            K=length(alphas);
            
            IS=zeros(length(alphas),1);
            
            for k=1:K
                
                alpha=alphas(k);
                
                w_k=alpha/2;
                
                % (1-lpha)x100 % prediction coverage
                
                Lt=localRowQuantile(curvesforecasts2(dataFperiod(i),:),alpha/2,weights2);
                Lt=(Lt>=0).*Lt;

                Ut=localRowQuantile(curvesforecasts2(dataFperiod(i),:),1-alpha/2,weights2);
                Ut=(Ut>=0).*Ut;

                coverage1=find(y>=Lt & y<=Ut);

                IS(k)=[(Ut-Lt)+(2/alpha).*(Lt-y).*(y<Lt)+(2/alpha)*(y-Ut).*(y>Ut)];

                sum1=sum1+w_k*IS(k);


            end

            %median prediction, m
            m=localRowQuantile(curvesforecasts2(dataFperiod(i),:),0.5,weights2);

            WISF(i)=(1/(K+1/2))*(w0*abs(y-m) + sum1);

           if isempty(weights2)
               quantilesf1(i,:)=quantile(curvesforecasts2(dataFperiod(i),:),alphaquantiles);
           else
               quantilesf1(i,:)=subepidemicWeightedQuantile(curvesforecasts2(dataFperiod(i),:),alphaquantiles,weights2);
           end
           quantilesf1(i,:)=(quantilesf1(i,:)>=0).*quantilesf1(i,:);


        end
    
        WISFS=[WISFS;[forecastingperiod1 mean(WISF)]];
        
    end
    
end

% -----------------------------------------------------------------------------
function q=localRowQuantile(curves,p,weights)
% Quantile p of the curves (columns) at each time point (row). Without weights
% this is the call made before the revision of 2026-11-24; with weights, the
% weighted quantile of the ensembles (subepidemicWeightedQuantile.m).
if isempty(weights)
    q=quantile(curves',p)';
else
    q=subepidemicWeightedQuantile(curves,p,weights);
end
