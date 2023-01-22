addpath('./TuckerL2E')
addpath('./TuckerL2E/L-BFGS-B-C/Matlab')
addpath('./TuckerL2E/tensor_toolbox-v3.2.1')
addpath('./RGrad')
addpath('./RGrad/functions')
%%
maxit = 100;

alpha0_set = 0.025:0.025:0.1;
alpha_set = [0.025,0.04,0.065,0.09];
% alpha_set = [0.05];
beta = 0.3;
gamma = 1;
d = [100,100,100];
m = length(d);
r = [2,2,2];
mu0 = 5;
sigma = 0.01;
% sigma_set = [0.01];
tol = 0.001;
repeat = 50;
S_amplitude = 0.1; % 0.1 small; 0.5 middle; 1 large

alg = 1; % 1 -- RGrad; 2 -- PGD
err = zeros(repeat,length(alpha_set));
err_l2e = zeros(repeat,length(alpha_set));
%% generate the low rank part
rng('default');
% Ttrue = randn(d);
% [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
% Ttrue = trimtensor(Ttrue,mu0); 

%% RGrad

fprintf('Using RGrad ......\n');
for k = 1:length(alpha_set)
    alpha0 = alpha0_set(k);
    alpha = alpha_set(k);
    for tt = 1:repeat
        fprintf('================================================\n');
        fprintf('Outlier sparsity = %f; repeating %d time(s)\n',alpha0,tt);

        %% generate random tensor 
        Ttrue = randn(d);
        [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
        %Ttrue = trimtensor(Ttrue,mu0);         
        Strue = S_amplitude * binornd(ones(d),alpha0) .* randn(d);
        Z = randn(d)*sigma;
        A = Ttrue + Strue + Z;
        fprintf('||T*||_F = %f, Frobenius norm of S* is: %f, maximum of S* is: %f, maximum of T* is: %f, ||Z||_F: %f\n',...
        norm(Ttrue(:)),norm(Strue(:)),max(abs(Strue(:))),max(abs(Ttrue(:))),norm(Z(:)));

        [T,~,time] = tRPCA_RGrad(A,r,gamma,alpha,Ttrue); % Ttrue only for calculating relative error

        err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        fprintf('Elapsed time is: %f; Relative error is: %f\n', time, norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
        [T,tau] = tucker_l2e_opt(tensor(A),2,'taumax',100,'init','HOSVD');
        err_l2e(tt,k) = norm(tensor(T)-tensor(Ttrue))/norm(tensor(Ttrue));
    end
end
mean(err,1)
mean(err_l2e,1)
csvwrite('err_rgrad.csv',err);
csvwrite('err_l2e.csv',err_l2e);


