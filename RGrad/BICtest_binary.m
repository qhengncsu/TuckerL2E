% testing BIC criterion for binary tensor 
clc;clear;
addpath(genpath(pwd));
%% set parameter
d = [100,100,100];
dprod = prod(d);
r0 = [3,3,3];
alpha0 = 0.01;
mu0 = 5;
gamma = 1;
zeta = 0.5;
stop_thres = 1e-2;
beta = 2;
kpr = 0.1;
maxit = 100;
sigma = 0.1;
repeat = 5;

%% generate A
rng('default');

Ttrue = randn(d);
[Ttrue,~,~,~,~,~] = hosvd(Ttrue,r0); 
Ttrue = trimtensor(Ttrue,mu0); Ttrue =  1*Ttrue/max(abs(Ttrue(:)));

Strue = binornd(ones(d),alpha0) .* randn(d); Strue = 0.1*Strue/max(abs(Strue(:)));


%%
p = @(x) 1./(1+exp(-x/sigma));
r_list = [1,2,3,4,5]; r_len = length(r_list);
alpha_list = 0.001:0.001:0.015; alpha_len = length(alpha_list);
% alpha_list = 0.002:0.002:0.02; alpha_list = [alpha_list(1:2) 0.005 alpha_list(3:end)]; alpha_len = length(alpha_list);
result = zeros(r_len,alpha_len,repeat);
fprintf('(alpha,r)              BIC              relative error         ||S||_0        ||S||_1          first term            -2log(L)\n')

for kk = 1:repeat
    A = binornd(ones(d),1./(1+exp(-(Ttrue + Strue)/sigma)));
    for i = 1:r_len
        r = r_list(i) * [1,1,1];
        for j = 1:alpha_len
            alpha = alpha_list(j);
            %% initialization
%             [T,C,U,UT] = init_binary(A,r,sigma);
            [T,C,U,UT,~,~] = hosvd(A,r); 
            S = threshold(T,gamma,alpha,A); 

            %% refinement
            Tprev = T;
            for k = 1:maxit
                G = -A./(1+exp((T+S)/sigma)) + (1-A)./(1+exp(-(T+S)/sigma));
                [G_p,~,~] = mani_proj(G,C,U,UT);
                W = T - beta*G_p;
                T = Trim2(W,r,zeta);
                S = gradPrune(T,gamma*alpha,kpr,A);

                if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres 
                    break;
                end
                Tprev = T;
            end
            result(i,j,kk) = bic_binary(A,S,T,r,sigma);
            fprintf('(%4f,%d)         %.4e        %.4e       %d          %.4e            %.6e           %.6e\n',alpha,r_list(i),...
               result(i,j),norm(T(:)-Ttrue(:))/norm(Ttrue(:)),nnz(S),sum(abs(S(:))),(nnz(S) + sum(d.*r))*log(dprod),-2*(sum(A.*log(p(T+S)),'all') + sum((1-A).*log(1-p(T+S)),'all')))
    %         fprintf('BIC result = %f, alpha = %f, r = %d, final error = %f\n',result(i,j),alpha,r_list(i),norm(T(:)-Ttrue(:))/norm(Ttrue(:)));
    %         fprintf('||S||_0 = %d,||A-T-S||_F = %f\n',nnz(S),norm(A(:)-T(:)-S(:)))
        end
    end
end
%% visualization
r1 = result(1,:);
r2 = result(2,:);
r3 = result(3,:);
r4 = result(4,:);
r5 = result(5,:);
%%
hold on
ax = gca;
ax.FontSize = 15; 
xticks(0.003:0.003:0.015);
xlabel('\alpha');
ylabel('BIC');
plot(alpha_list,r1,'r-s');
plot(alpha_list,r2,'g-o');
plot(alpha_list,r3,'b-*','LineWidth',4);
plot(alpha_list,r4,'m-+');
plot(alpha_list,r5,'k-^');


grid minor
% legend({'r = 3'});
legend({'r = 1','r = 2','r = 3','r = 4','r = 5'});
legend('Location','Best');
hold off
    
    
    