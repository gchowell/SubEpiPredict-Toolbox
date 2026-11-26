
function [LB,UB]=getbounds(npatches,data)

global flag1

r_bnds=[];
p_bnds=[];
a_bnds=[];
K_bnds=[];

rlb=mean(abs(data(1:2,1)))/200;
rub=max(abs(data(1:2,1)))*50;

%rlb
%rub

for j=1:npatches
    
    r_bnds=[r_bnds [rlb;rub]];
    
    Kmax=10000000000;
    
    switch flag1
        
        case 0 %GGM
            
            p_bnds=[p_bnds [0;1]];
            a_bnds=[a_bnds [1;1]];
            
            K_bnds=[K_bnds [20;20]];
            
            
        case 1 %GLM
            
            p_bnds=[p_bnds [0;1]];
            a_bnds=[a_bnds [1;1]];
            K_bnds=[K_bnds [20;Kmax]];
            
        case 2 %GRM
            p_bnds=[p_bnds [0;1]];
            a_bnds=[a_bnds [0;10]];
            K_bnds=[K_bnds [20;Kmax]];
            
            
        case 3 %Logistic
            p_bnds=[p_bnds [1;1]];
            a_bnds=[a_bnds [1;1]];
            K_bnds=[K_bnds [20;Kmax]];
            
        case 4 %Richards
            p_bnds=[p_bnds [1;1]];
            a_bnds=[a_bnds [0;10]];
            K_bnds=[K_bnds [20;Kmax]];
            
        case 5 %Gompertz
            p_bnds=[p_bnds [1;1]];
            a_bnds=[a_bnds [0;10]];
            K_bnds=[K_bnds [20;Kmax]];
            
    end
end



LB=[r_bnds(1,:) p_bnds(1,:) a_bnds(1,:) K_bnds(1,:)]; %r p a K alpha d

UB=[r_bnds(2,:) p_bnds(2,:) a_bnds(2,:) K_bnds(2,:)]; %r p a K alpha d

