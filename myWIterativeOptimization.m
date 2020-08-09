function [Wnew,Obj] = myWIterativeOptimization(W,K0,T,mis_indx)
%% Kp: n*n;
%% T: n*n: 
%% obs_indx: index set of observed samples
%% W: c*m;
n = size(T,1);
obs_indx = setdiff(1:n,mis_indx);
Tcm = T(obs_indx,mis_indx);
Tmm = T(mis_indx,mis_indx);
flag = 1;
iter = 0;
while flag
    iter = iter +1;
%     WK = W'*K0;%% m*c
    %% WK' = K0*W;
%     Grad1 = 2*K0*WK';
%     Grad2 = 2*K0*Tcm;
%     Grad3 = 2*(WK'*W')*WK';
%     Grad4 = 2*WK'*Tmm;
%     Grad13 = Grad1+Grad3;
%     Grad24 = Grad2+Grad4;
    Wnew = ((eye(length(obs_indx))+W*W')*K0)\(Tcm+W*Tmm);
    
    WnK = Wnew'*K0;
    Obj1 = trace(WnK*WnK');
    Obj2 = trace(WnK*Tcm);
    Obj3 = trace((WnK*Wnew)*(WnK*Wnew));
    Obj4 = trace((WnK*Wnew)*Tmm);
    Obj(iter) = Obj1 - 2*Obj2 + (1/2)*Obj3 - Obj4;
    if (iter>2 && (Obj(iter-1)-Obj(iter))/Obj(iter)<1e-3)
        flag =0;
    end
    W = Wnew;    
end