clc;clear;
addpath(genpath('functions'));
%% set parameters
param.d = [100,100,100];
param.r = [2,2,2];

param.alpha = 0.02;
param.gamma = 2;
param.beta = 0.5; % stepsize
param.mu0 = 5; % incoherent parameter
param.sigma = 0.01; % noise level
param.perturb_sigma = 0.1;
param.init_hosvd = 1;
param.stop_thres = 0.001;

%% modify alpha
% alpha_list = [0.001,0.005,0.01,0.05,0.1];
% beta_list = [0.25,0.25,0.35,0.5,0.5];
% result = cell(length(alpha_list),1);
% for i = 1:length(alpha_list)
%     param.beta = beta_list(i);
%     param.alpha = alpha_list(i)
%     result{i} = tensorRPCA(param);
% end
% save(['result/tRPCA/test_modify_alpha/result',datestr(datetime),'.mat'],'result')

%% modify sigma
sigma_list = [0.001,0.005,0.01,0.02,0.03];
stop_thers_list = [0.0001,0.0001,0.0001,0.005,0.005];
beta_list = [0.2,0.2,0.2,0.2,0.6];
% sigma_list = [0.03];
% stop_thers_list = [0.005];
% beta_list = [0.3];

result = cell(length(sigma_list),1);
for i = 1:length(sigma_list)
    param.sigma = sigma_list(i);
    param.stop_thres = stop_thers_list(i);
    param.beta = beta_list(i)
    result{i} = tensorRPCA(param);
end
save(['result/tRPCA/test_modify_sigma/result',datestr(datetime),'.mat'],'result')

%% modify gamma
gamma_list = [2,5,10];
stop_thers_list = [0.001,0.001,0.001];
beta_list = [0.5,0.5,0.5];
result = cell(length(gamma_list),1);
for i = 1:length(gamma_list)
    param.gamma = gamma_list(i);
    param.stop_thres = stop_thers_list(i);
    param.beta = beta_list(i);
    result{i} = tensorRPCA(param);
end
% save(['result/tRPCA/test_modify_gamma/result',datestr(datetime),'.mat'],'result')

%%
r1 = result{1}.relerrList;
r2 = result{2}.relerrList;
r3 = result{3}.relerrList;

hold on;

plot(log(r1(13:end)),'r-o');
plot(log(r2(13:end)),'b-+');
plot(log(r3(13:end)),'m-^');

ax = gca;
ax.FontSize = 15; 


xlabel('iterations')
ylabel('log(relative error)')
grid minor;
legend({'\gamma = 2','\gamma = 5','\gamma = 10'});
legend('Location','Best');
hold off;