function objfunction=plotModifiedLogisticGrowthPatchMethodsFullLogLik(z)

% Use the same initial-value problem for calibration and reconstruction.
global flag1 method1 timevect ydata I0 npatches onset_thr yfit
global onset_fixed

% Reject stale or mis-sized parameter vectors BEFORE any block indexing.
[rs1,ps1,as1,Ks1,alpha,d]=unpackSubepidemicParameters(z,npatches);

[~,~,yfit]=simulateSubepidemic(timevect,I0,npatches,onset_fixed, ...
    onset_thr,flag1,rs1,ps1,as1,Ks1);

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

        case 6 % Sum of Absolute Deviations (SAD)

            objfunction=sum(abs(ydata-yfit));


    end

end


