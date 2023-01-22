clc;clear;
addpath(pwd);
%%
alpha = 0.0001;
gamma = 1.5;
d = [100,100,100];
r = [2,2,2];
mu0 = 5;
sigma = 0.00;

%% generate the observed tensor A
rng('default');
Ttrue = randn(d);
[Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
Ttrue = trimtensor(Ttrue,mu0); 

Strue = binornd(ones(d),alpha) .* randn(d);
Z = randn(d)*sigma;
A = Ttrue + Strue + Z;
unfold_A = ndim_unfold(A,1);


%%
params.max_iter = 100;
params.tol = 0.001;
params.step_const = 0.5;

tic;
[U,V] = rpca_gd(unfold_A, 2, alpha*gamma, params);
T = U*V';
T = ndim_fold(T,1,d);
norm(T(:)-Ttrue(:))/norm(Ttrue(:))
toc;