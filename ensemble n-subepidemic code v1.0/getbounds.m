
function [LB,UB]=getbounds(npatches)

global flag1

r_bnds=[];
p_bnds=[];
a_bnds=[];
K_bnds=[];

for j=1:npatches
    
    r_bnds=[r_bnds [0;2000]];
    
    Kmax=100000000000;
    
    switch flag1(j)
        
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

