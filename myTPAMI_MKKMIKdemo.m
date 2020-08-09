clear
clc
warning off;

path = 'D:\study\code\work2016\';
addpath(genpath(path));
dataName = 'bbcsport2view'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
%% %% washington; wisconsin; texas; cornell
load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;
% [H_normalized,gamma,obj] = mkkmeans_train(KH,numclass,qnorm);
% res_gnd = myNMIACC(H_normalized,Y,numclass);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsionset = [0.1:0.1:0.9];
for ie = 1:length(epsionset)
     for iter = 1 
        load([path,'generateAbsentMatrix\',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_missingIndex_iter_',num2str(iter),'.mat'],'S');
        mis_indx = S{1}.indx; 
        obs_indx = setdiff(1:num,mis_indx);
        num_mis = length(mis_indx);
        num_obs = num - num_mis;
        W = rand(num_obs,num_mis);
        Kp = KH(obs_indx,obs_indx,1);
        T = rand(num);
        T = (T+T')/2;
        [Wnew,Obj] = myWIterativeOptimization(W,Kp,T,mis_indx);
      end
end