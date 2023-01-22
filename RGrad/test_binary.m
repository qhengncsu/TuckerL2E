clc;clear;
addpath(genpath('functions'));
%%
param.d = [100,100,100];
param.r = [2,2,2];

%%
param.mu0 = 5;
param.beta = 2; %
param.sigma = 1; %
param.zeta = 18;
param.perturb_sigma = 10;
param.init_hosvd = 0;
param.gamma = 1.1;
param.alpha = 0.001; 
param.kpr = 10;
param.scaling_param_T = 200;
param.scaling_param_S = 3;
param.stop_thres = 0.005;
param.randseed = 1;

default_param = param;
%% modify alpha
param = default_param;
alpha_list = [0.001,0.0050,0.01,0.02];
result = cell(length(alpha_list),1);
for i = 1:length(alpha_list)
    param.alpha = alpha_list(i);
    result{i} = binary(param);
end
% save(['result/binary/modify_alpha/result',datestr(datetime),'.mat'],'result');

%% 
r1 = result{1}.relerrList;
r2 = result{2}.relerrList;
r3 = result{3}.relerrList;
r4 = result{4}.relerrList;

hold on;
plot(log(r1),'m-^');
plot(log(r2),'k-+');
plot(log(r3),'b-o');
% plot(r4,'r-s');

legend({'\alpha = 0.001','\alpha = 0.005','\alpha = 0.01'})
legend('Location','Best')
xlabel('iterations')
ylabel('log(relative error)')
grid minor;
set(gca,'Fontsize',15);
hold off;

%% modify sigma
% param = default_param;
% sigma_list = [0.8,0.9,1,1.1,1.2];
% result = cell(length(sigma_list),1);
% for i = 1: length(sigma_list)
%     param.sigma = sigma_list(i)
%     result{i} = binary(param);
% end
% save(['result/binary/modify_sigma/result',datestr(datetime),'.mat'],'result');
%% modify gamma
% param = default_param;
% gamma_list = [1,2,3,4,5,6,7,8];
% result = cell(length(gamma_list),1);
% for i = 1:length(gamma_list)
%     param.gamma = gamma_list(i)
%     result{i,1} = binary(param);
% end
% save(['result/binary/modify_gamma/result',datestr(datetime),'.mat'],'result');
%% check for the best param. when alpha = 0.02
% param = default_param;
% param.alpha = 0.02;
% beta_list = [0.05,0.1,0.25,0.5,0.8];
% for j = 1:length(beta_list)
%     param.beta = beta_list(j)
%     result = binary(param);
% end

%% check for the best param. when alpha = 0.03
% param = default_param;
% param.alpha = 0.03;
% beta_list = [0.05,0.1,0.25,0.5,0.8];
% for j = 1:length(beta_list)
%     param.beta = beta_list(j)
%     result = binary(param);
% end

%% error bar -- fix sigma
% gamma_list = [1,1.5,2,2.5,3,3.5,4];
% alpha_list = [0.001,0.002,0.005,0.01,0.015,0.02];
% % gamma_list = [1,2];
% % alpha_list = [0.001];
% for i = 1: length(gamma_list)
%     param.gamma = gamma_list(i);
%     result = cell(length(alpha_list),5);
%     for j = 1:length(alpha_list)
%         param.alpha = alpha_list(j);
%         for k =1 :5
%             param.randseed = k*10 + 2
%             result{j,k} = binary(param);
%         end
%     end
%     save_dir = ['result/binary/errBar_fix_sigma/resultGamma = ',num2str(param.gamma),'--',datestr(datetime),'.mat'];
%     save(save_dir,'result');
% end

%% error bar -- fix gamma
sigma_list = [0.8,0.9,0.95,1,1.1,1.2];
alpha_list = [0.001,0.002,0.005,0.01,0.015,0.02];
% sigma_list = [0.8,0.9];
% alpha_list = [0.001];
for i = 1: length(sigma_list)
    param.sigma = sigma_list(i);
    result = cell(length(alpha_list),5);
    for j = 1:length(alpha_list)
        param.alpha = alpha_list(j);
        for k =1 :5
            param.randseed = k*10 + 2
            result{j,k} = binary(param);
        end
    end
    save_dir = ['result/binary/errBar_fix_gamma/resultSigma = ',num2str(param.sigma),'--',datestr(datetime),'.mat'];
    save(save_dir,'result');
end



