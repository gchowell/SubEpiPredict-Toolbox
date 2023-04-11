function objfunction=plotModifiedLogisticGrowthPatchMethods1(z)

global flag1 method1 timevect ydata I0 npatches onset_thr yfit

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed


invasions=zeros(npatches,1);
timeinvasions=zeros(npatches,1);
Cinvasions=zeros(npatches,1);

invasions(1)=1;
timeinvasions(1)=0;
Cinvasions(1)=0;

rs1=z(1:npatches);
ps1=z(npatches+1:2*npatches);
as1=z(2*npatches+1:3*npatches);
Ks1=z(3*npatches+1:4*npatches);

alpha=z(end-1);
d=z(end);

%if onset_thr>Ks1(1)

%    npatches=1;
%end


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


[t,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs1,ps1,as1,Ks1,npatches,onset_thr,flag1);

% if sum(invasions)<npatches
%
%     npatches=sum(invasions);
%
%     z=[rs_hat(1:npatches) ps_hat(1:npatches) as_hat(1:npatches) Ks_hat(1:npatches) alpha_hat d_hat];
%
%     'entro 1'
%
%     pause
%
% end

y=sum(x,2);

totinc=[y(1,1);diff(y(:,1))];

if onset_fixed==0
    totinc(1)=totinc(1)-(npatches-1);
end

yfit=totinc;


eps=0.001;

%%MLE expression
%This is the negative log likelihood, name is legacy from least squares code
%Note that a term that is not a function of the params has been excluded so to get the actual
%negative log-likliehood value you would add: sum(log(factorial(sum(casedata,2))))


if sum(yfit)==0
    objfunction=10^10;%inf;
else
    %    z
    yfit(yfit==0)=eps; %set zeros to eps to allow calculation below.  Shouldn't affect solution, just keep algorithm going.
    
    switch method1
        
        case 0

            %minimize sum of square errors
            objfunction=sum((ydata-yfit).^2);
            
            %without constants
            %objfunction= nd*log(SSE);   % MLE normal distribution (Least squares)
            
            %full negative log-likelihood:
            %objfunction= -1*((-length(ydata)/2)*log(2*pi)-(length(ydata)/2)*log(SSE/length(ydata))-length(ydata)/2);
            
            
        case 1 %MLE Poisson (negative log-likelihood)
            
            
            %             sum1=0;
            %             for i=1:length(ydata)
            %
            %                 sum1=sum1+ydata(i)*log(yfit(i))-sum(log(2:1:ydata(i)))-yfit(i);
            %
            %             end
            %
            %             objfunction=-sum1;
            
            %without constants
            objfunction=-sum(ydata.*log(yfit)-yfit);
            
            
        case 3  % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;
            
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha)*yfit(i));
                    
                end
                
                %sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha)-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha);

            end
            
            objfunction=-sum1;
            
            
        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha));
                    
                end
                
                %sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i))-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i));

            end
            
            objfunction=-sum1;
            
            
        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha)*yfit(i).^(2-d));
                    
                end
                
                %sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-1)))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-1)))-sum(log(2:1:ydata(i)));
                sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-2)).*yfit(i))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-2)).*yfit(i));

            end
            
            objfunction=-sum1;
            
    end
    
end


