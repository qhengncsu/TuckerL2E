% For different values of alpha0, we use BIC to estimate alpha0
% then we repeat several times using the estimated alpha

clc;clear;
addpath(genpath('functions'));
%%
maxit = 100;
alpha0 = 0.1;
alpha_list = 0.075:0.005:0.125; alpha_len = length(alpha_list);
beta = 0.3;
gamma = 1;
d = [100,100,100];
m = length(d);
r = [2,2,2];
mu0 = 5;
sigma = 0.01;
tol = 0.001;
S_amplitude = 0.5; % 0.1 small; 0.5 middle; 1 large

%% generate the low rank part
rng('default');

%% generate random tensor 
Ttrue = randn(d);
[Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
Ttrue = trimtensor(Ttrue,mu0);         
Strue = S_amplitude * binornd(ones(d),alpha0) .* randn(d);
Z = randn(d)*sigma;
A = Ttrue + Strue + Z;
for k = 1:alpha_len
    alpha = alpha_list(k);
    
%     fprintf('||T*||_F = %f, Frobenius norm of S* is: %f, maximum of S* is: %f, maximum of T* is: %f, ||Z||_F: %f\n',...
%     norm(Ttrue(:)),norm(Strue(:)),max(abs(Strue(:))),max(abs(Ttrue(:))),norm(Z(:)));

    [T,S,time] = tRPCA_RGrad(A,r,gamma,alpha,Ttrue); % Ttrue only for calculating relative error

    fprintf('alpha = %.4f; BIC value: %.6e\n', alpha, bic(A,S,T,r));
end





