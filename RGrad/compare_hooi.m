clc;clear;
addpath(genpath('functions'));
%%
maxit = 100;

alpha_set = 0.025:0.025:0.1;
% alpha_set = [0.025];
d = [100,100,100];
m = length(d);
r = [2,2,2];
mu0 = 5;
sigma = 0.01;
S_amplitude = 1; 
repeat = 7;

err = zeros(repeat,length(alpha_set));

%%
rng('default');

%%
for k = 1:length(alpha_set)
    alpha = alpha_set(k);
    for tt = 1:repeat
        fprintf('Current sparsity: %f, repeating %d times\n',alpha,tt);
        %% generate random tensor 
        Ttrue = randn(d);
        [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
        Ttrue = trimtensor(Ttrue,mu0); 
        Strue = S_amplitude * binornd(ones(d),alpha) .* randn(d);
        Z = randn(d)*sigma;
        A = Ttrue + Strue + Z;
        
        tic;
        [T] = HOOI(A,r);
        time = toc;
        
        err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        fprintf('Elapsed time is: %f; Relative error is: %f\n', time, norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
    end
end

%%
result.err = err;
result.S_amplitude = S_amplitude;
result.sigma = sigma;
result.alpha_set = alpha_set;
result.alg = 5;

