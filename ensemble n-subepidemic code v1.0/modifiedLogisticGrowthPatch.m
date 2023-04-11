function dx=modifiedLogisticGrowthPatch(t,x,rs1,ps1,as1,Ks1,npatches,onset_thr,flag1)

global invasions
global timeinvasions
global Cinvasions

dx=zeros(npatches,1);

for j=1:npatches
    
    if invasions(j)==0
        
        %K1=K*exp(-(j-1).^q);
        
        %if K1<onset_thr
        %    K1=onset_thr+0.1;
        %end
        
        invasions(j)=(x(j-1,1)>=onset_thr);
        timeinvasions(j)=invasions(j).*t;
        Cinvasions(j)=invasions(j).*x(j-1,1);
        
    end
    
    switch flag1
        
        case 0 %GGM
            
            dx(j,1)=invasions(j).*rs1(j)*x(j,1).^ps1(j);
            
        case 1 % GLM
            
            dx(j,1)=invasions(j).*(rs1(j)*(x(j,1).^ps1(j))*(1-(x(j,1)/Ks1(j))));
            
        case 2 % Gen Richards
            dx(j,1)=invasions(j).*(rs1(j)*(x(j,1).^ps1(j))*(1-(x(j,1)/Ks1(j))).^as1(j));
            
        case 3 %Logistic
            dx(j,1)=invasions(j).*(rs1(j)*(x(j,1))*(1-(x(j,1)/Ks1(j))));
            
        case 4 %Richards
            dx(j,1)=invasions(j).*(rs1(j)*(x(j,1))*(1-(x(j,1)/Ks1(j)).^as1(j)));
            
        case 5 %Gompertz
            %dx(j,1)=invasions(j).*(rs1(j)*x(j,1)*exp(-as1(j)*t));
            dx(j,1)=invasions(j).*(rs1(j)*x(j,1).*(log(Ks1(j)/x(j,1))));
            %dx(j,1)=invasions(j).*(rs1(j)*x(j,1).*(log(Ks1(j)/x(j,1))).^as1(j));
            
    end
    
end
