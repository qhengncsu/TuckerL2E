clc;clear;
addpath(genpath('functions'));
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
repeat = 5;
S_amplitude = 0.1; % 0.1 small; 0.5 middle; 1 large

alg = 1; % 1 -- RGrad; 2 -- PGD
err = zeros(repeat,length(alpha_set));
%% generate the low rank part
rng('default');
% Ttrue = randn(d);
% [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
% Ttrue = trimtensor(Ttrue,mu0); 

%% RGrad
if alg == 1
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
            Ttrue = trimtensor(Ttrue,mu0);         
            Strue = S_amplitude * binornd(ones(d),alpha0) .* randn(d);
            Z = randn(d)*sigma;
            A = Ttrue + Strue + Z;
            fprintf('||T*||_F = %f, Frobenius norm of S* is: %f, maximum of S* is: %f, maximum of T* is: %f, ||Z||_F: %f\n',...
            norm(Ttrue(:)),norm(Strue(:)),max(abs(Strue(:))),max(abs(Ttrue(:))),norm(Z(:)));
        
            [T,~,time] = tRPCA_RGrad(A,r,gamma,alpha,Ttrue); % Ttrue only for calculating relative error
            
            err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
            fprintf('Elapsed time is: %f; Relative error is: %f\n', time, norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
        end
    end
    result.err = err;
    result.alg = alg;
    result.S_amplitude = S_amplitude;
    result.alpha_set = alpha_set;
    result.sigma = sigma;
    
%% PGD
elseif alg == 2
    fprintf('Using PGD ......\n');
    for k = 1:length(alpha_set)
        alpha = alpha_set(k);
        for tt = 1:repeat
            fprintf('================================================\n');
            fprintf('Outlier sparsity = %f; repeating %d time(s)\n',alpha,tt);
            %% generate random tensor 
            Ttrue = randn(d);
            [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
            Ttrue = trimtensor(Ttrue,mu0); 
            
            Strue = S_amplitude * binornd(ones(d),alpha) .* randn(d);

            Z = randn(d)*sigma;
            A = Ttrue + Strue + Z;

            %% initialization 
            [T,C,U,UT,~,~] = hosvd(A,r); % use HOSVD as initialization
            init_err = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
            disp(['initialization rel err.: ', num2str(init_err)])
            tic;
            for i = 1:maxit
                Tprev = T;
                G = T - A;
                W = T - beta*G;
                [T,C,U,UT,~,~] = hosvd(W,r); 
                T = trimtensor(T,mu0);
                disp(['iter:', num2str(i),', rel err.:',num2str(norm(T(:) - Ttrue(:))/norm(Ttrue(:)))])
                if norm(T(:) - Tprev(:))/norm(T(:)) < 0.001
                    break;
                end
            end
            time = toc;
            err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
            fprintf('Elapsed time is: %f; Relative error is: %f\n', time, norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
        end
    end
    result.err = err;
    result.alg = alg;
    result.S_amplitude = S_amplitude;
    result.alpha_set = alpha_set;
    result.sigma = sigma;
end



