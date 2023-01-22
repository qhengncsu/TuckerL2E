% This file is the comparison experiment between RGrad and PGD


clc;clear;
addpath(genpath('functions'));

%% parameter setting - fixed parameters
maxit = 100;
d = [100,100,100];
m = length(d);
r = [2,2,2];
tol = 0.01;

% tuning parameter
beta = 3;
gamma = 1;
alpha0_set = 0.005:0.005:0.02;
alpha_set = [0.003,0.007,0.012,0.016];
% alpha_set = [0.025];
zeta = 2.5;
S_amp = 10;
T_amp = 5;
kpr = 1;
mu0 = 5;
sigma = 5;

% other parameters
repeat = 5;
err = zeros(repeat,length(alpha_set));
alg = 1;

%% model
f       = @(x) (1 ./ (1 + exp(-x/sigma)));
fprime  = @(x) (1/sigma) * (exp(x/sigma) ./ (1 + exp(x/sigma)).^2);  

%%
rng('default');
if alg == 1
%     Ttrue = T_amp * randn(d);
%     [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
%     Ttrue = trimtensor(Ttrue,mu0); 
    Ttrue = T_amp*generate_low_rank(d,r);
    for k = 1:length(alpha_set)
        alpha0 = alpha0_set(k);
        alpha = alpha_set(k);
        for tt = 1:repeat
            fprintf('================================================\n');
            fprintf('Outlier sparsity = %f; repeating %d time(s)\n',alpha0,tt);
            %% generate the observed tensor
            
%             Strue = S_amp * binornd(ones(d),alpha) .* randn(d);
            Strue = S_amp * binornd(ones(d),alpha0);
            A = binornd(ones(d),f(Strue + Ttrue));

            %% initialization
%             [T,C,U,UT,~,~] = hosvd(A,r); 
            [T,C,U,UT] = init_binary(A,r,sigma);
            S = gradPrune(T,gamma*alpha,kpr,A);
            fprintf('initialization rel err.: %f\n', norm(T(:) - Ttrue(:))/norm(Ttrue(:)))

            %% iteration
            for i = 1:maxit
                Tprev = T;
                G = -A.*fprime(T+S)./f(T+S) + (1-A).*fprime(T+S)./(1-f(T+S));
                [G,~,~] = mani_proj(G,C,U,UT);
                W = T - beta*G;
                [T,C,U,UT] = Trim2(W,r,zeta);
                S = gradPrune(T,gamma*alpha,kpr,A);
                fprintf('relative error is %f\n',norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
                if norm(T(:) - Tprev(:))/norm(T(:)) < tol
                    break;
                end
            end
            err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        end
    end
    
%% PGD
elseif alg == 2
%     Ttrue = T_amp * randn(d);
%     [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
%     Ttrue = trimtensor(Ttrue,mu0); 
    Ttrue = T_amp*generate_low_rank(d,r);
    for k = 1:length(alpha_set)
        alpha = alpha_set(k);
        for tt = 1:repeat
            fprintf('================================================\n');
            fprintf('Outlier sparsity = %f; repeating %d time(s)\n',alpha,tt);
            %% generate the observed tensor
%             Strue = S_amp * binornd(ones(d),alpha) .* randn(d);
            Strue = S_amp * binornd(ones(d),alpha);
            A = binornd(ones(d),f(Strue + Ttrue));

            %% initialization
%             [T,C,U,UT,~,~] = hosvd(A,r); 
            [T,C,U,UT] = init_binary(A,r,sigma); 
            fprintf('initialization rel err.: %f\n', norm(T(:) - Ttrue(:))/norm(Ttrue(:)))

            %% iteration
            for i = 1:maxit
                Tprev = T;
                G = -A.*fprime(T)./f(T) + (1-A).*fprime(T)./(1-f(T));
                W = T - beta*G;
                [T,C,U,UT] = hosvd(W,r);
                fprintf('relative error is %f\n',norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
                if norm(T(:) - Tprev(:))/norm(T(:)) < tol
                    break;
                end
            end
            err(tt,k) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        end
    end
end
