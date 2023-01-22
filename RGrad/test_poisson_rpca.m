clc;clear;
addpath(genpath(pwd));
%% parameters
d = [100,100,100]; dprod = prod(d);
r = [2,2,2];
alpha = 0.001;
mu0 = 5;
gamma = 1.1;
stop_thres = 1e-3;
beta = 0.1;

zeta = 0.5;
kpr = 0.5;

maxit = 100;
I_list = [10,20,50,100]; Ilen = length(I_list);
repeat = 5;
%% generate A
rng('default');
Ttrue = randn(d); [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
Ttrue = trimtensor(Ttrue,mu0); 
Ttrue =  0.5*Ttrue/max(abs(Ttrue(:)));

Strue = binornd(ones(d),alpha) .* randn(d); Strue = 0.5*Strue/max(abs(Strue(:)));

%% 
err = zeros(Ilen,repeat);


%% refinement
for i = 1:Ilen
    I = I_list(i);
    for j = 1:repeat
        A = poissrnd(I*exp(Ttrue + Strue));
        
        %% initialization
        tildeT = log((A + 0.5)/I);
        [T,C,U,UT,~,~] = hosvd(tildeT,r);
        fprintf('initialization error: %.4f\n', norm(T(:)-Ttrue(:))/norm(Ttrue(:)));
        S = gradPrune_pois(T,gamma*alpha,kpr,A,I);

        Tprev = T;
        for k = 1:maxit
            G = -A/I + exp(T+S);
            [G_p,~,~] = mani_proj(G,C,U,UT);
            W = T - beta*G_p;
            T = Trim2(W,r,zeta);
            S = gradPrune_pois(T,gamma*alpha,kpr,A,I);

%             fprintf('step %d, relative error %.4f\n',k,norm(T(:)-Ttrue(:))/norm(Ttrue(:)));

            if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres 
                break;
            end
            Tprev = T;
        end
        fprintf('final err: %.4f\n',norm(T(:)-Ttrue(:))/norm(Ttrue(:)));
        err(i,j) = norm(T(:)-Ttrue(:))/norm(Ttrue(:));
    end
end





