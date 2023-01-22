% For different values of alpha0, we use BIC to estimate alpha0
% then we repeat several times using the estimated alpha

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
alpha0 = 0.01;
alpha_list = 0.001:0.001:0.03; alpha_len = length(alpha_list);
zeta = 2.5;
S_amp = 1;
T_amp = 5;
kpr = 1;
mu0 = 5;
sigma = 5;

%% model
f       = @(x) (1 ./ (1 + exp(-x/sigma)));
fprime  = @(x) (1/sigma) * (exp(x/sigma) ./ (1 + exp(x/sigma)).^2);  

%% 
rng('default');
%     Ttrue = T_amp * randn(d);
%     [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
%     Ttrue = trimtensor(Ttrue,mu0); 
Ttrue = T_amp*generate_low_rank(d,r);
%% generate the observed tensor
Strue = S_amp * binornd(ones(d),alpha0);
A = binornd(ones(d),f(Strue + Ttrue));



%% iteration
for k = 1:alpha_len
    
    %% initialization
    alpha = alpha_list(k);

    [T,C,U,UT,~,~] = hosvd(A,r); 
    % [T,C,U,UT] = init_binary(A,r,sigma);
    S = gradPrune(T,gamma*alpha,kpr,A);
    % fprintf('initialization rel err.: %f\n', norm(T(:) - Ttrue(:))/norm(Ttrue(:)))
    
    for i = 1:maxit
        Tprev = T;
        G = -A.*fprime(T+S)./f(T+S) + (1-A).*fprime(T+S)./(1-f(T+S));
        [G,~,~] = mani_proj(G,C,U,UT);
        W = T - beta*G;
        [T,C,U,UT] = Trim2(W,r,zeta);
        S = gradPrune(T,gamma*alpha,kpr,A);
    %     fprintf('relative error is %f\n',norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
        if norm(T(:) - Tprev(:))/norm(T(:)) < tol
            break;
        end
    end

    fprintf('alpha = %.4f; BIC value: %.6e\n', alpha, bic_binary(A,S,T,r,sigma));
end 
