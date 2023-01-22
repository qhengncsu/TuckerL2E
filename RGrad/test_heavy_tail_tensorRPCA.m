clc;clear;
addpath(genpath('functions'));
%% set parameters
param.d = [100,100,100];
param.r = [2,2,2];

% params. for generating Z
param.theta = 5;
param.scaling_z = 0.01;
param.sigma_z = sqrt(param.theta/(param.theta-2))*param.scaling_z;
% params. for initialization
param.perturb_sigma = 0.05;
param.init_hosvd = 0;
% params. for iteration
param.beta = 0.5; 
param.stop_thres = 0.001;
% params. for updating S
param.alpha = 0.01;
param.gamma = 2;
param.mu0 = 5; 

%% modify theta
% theta_list = [2.2,2.5,3,3.5];
% stop_thers_list = [0.001,0.001,0.001,0.001];
% alpha_list = [.01,.01,.01,.01];
% beta_list = [0.5,0.5,0.5,0.5];
% 
% result = cell(length(theta_list),1);
% for i = 1:length(theta_list)
%     param.theta = theta_list(i);    
%     param.stop_thres = stop_thers_list(i);
%     param.alpha = alpha_list(i);
%     param.beta = beta_list(i)
%     result{i} = heavy_tail_tensorRPCA(param);
% end
% save(['result/heavy_tail_tRPCA/test_modify_theta/result',datestr(datetime),'.mat'],'result')

%% modify alpha
% alpha_list = [0.0,0.01,0.02,0.03,0.04];
% stop_thers_list = [0.000,0.0001,0.0001,0.0001,0.0001];
% beta_list = [.7,.7,.7,.7,.7];
% % alpha_list = [0.01];
% % stop_thers_list = [0.0];
% % beta_list = [0.7];
% result = cell(length(alpha_list),1);
% for i = 1:length(alpha_list)
%     param.alpha = alpha_list(i);
%     
%     param.stop_thres = stop_thers_list(i);
%     param.beta = beta_list(i)
%     result{i} = heavy_tail_tensorRPCA(param);
% end
% save(['result/heavy_tail_tRPCA/test_modify_alpha/result',datestr(datetime),'.mat'],'result')

%% errBar: modify alpha
theta_list = [2.05,2.1,2.2,2.5,3,4,5];
alpha_list = [0.01,0.02,0.05,0.1,0.15,0.2];
stop_thers_list = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001];
beta_list = [.5,.5,.5,.5,.5,.5];
% theta_list = [2.05,2.1]; 
% alpha_list = [0.0];
% stop_thers_list = [0.0001];
% beta_list = [0.5];
result = cell(length(theta_list),length(alpha_list),5);
for l = 1:length(theta_list)
    param.theta = theta_list(l);
    result = cell(length(alpha_list),5);
    for i = 1:length(alpha_list)
        param.alpha = alpha_list(i);
        param.sigma_z = sqrt(param.theta/(param.theta-2))*param.scaling_z;

        param.stop_thres = stop_thers_list(i);
        param.beta = beta_list(i);
        for j = 1:5
            param.randseed = j*10 + 2
            result{i,j} = heavy_tail_tensorRPCA(param);
        end
    end
    save_dir = ['result/heavy_tail_tRPCA/errBar_modify_alpha/resultTheta = ',num2str(param.theta),'--',datestr(datetime),'.mat'];
%     save(save_dir,'result')
end




