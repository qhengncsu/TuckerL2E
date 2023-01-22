% testing BIC criterion
clc;clear;
addpath(genpath(pwd));
%% set parameter
d = [100,100,100];
dprod = prod(d);
r0 = [3,3,3];
alpha0 = 0.05;
mu0 = 5;
sigma = 5e-3;
gamma = 1;
stop_thres = 1e-3;
beta = 0.3;
maxit = 100;

%% generate A
load 'data/it4tnsr_notime.mat';
A = log(1+it4tnsr_notime);

rng('default');
per = 0.1;
[A,testA,nnz_idx] = data_partition(A,per);

%%
r_list = [1,2,3]; r_len = length(r_list);
alpha_list = 0.002:0.002:0.07; alpha_len = length(alpha_list);
result = zeros(r_len,alpha_len);
fprintf('(alpha,r)              BIC              relative error             first term             second term\n')

for i = 1:r_len
    r = r_list(i) * [1,1,1];
    for j = 1:alpha_len
        alpha = alpha_list(j);
        %% initialization
        [T,C,U,UT,~,~] = hosvd(A,r); 
        S = threshold(T,gamma,alpha,A); 
%         fprintf('initialization error: %f\n', norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
        %% refinement
        Tprev = T;
        for k = 1:maxit
            G = T + S - A;
            [G_p,~,~] = mani_proj(G,C,U,UT);
            W = T - beta*G_p;
            [T,C,U,UT,~,~] = hosvd(W,r); 
            tildeT = trimtensor(T,mu0); 
            S = threshold(tildeT,gamma,alpha,A); 

            if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres
                break;
            end
            Tprev = tildeT;
        end
        result(i,j) = bic(A,S,T,r);
        err = (T+S).*nnz_idx - testA;
        fprintf('(%4f,%d)         %.4e        %.4e                 %.4e           %.4e\n',alpha,r_list(i),...
                       result(i,j),norm(err(:)),(nnz(S) + sum(d.*r))*log(dprod),dprod*log(norm(A(:)-T(:)-S(:))^2))
    end
end

%% visualization
% r1 = log(result(1,:));
% r2 = log(result(2,:));
% r3 = log(result(3,:));
r1 = result(1,:);
r2 = result(2,:);
r3 = result(3,:);

hold on
ax = gca;
ax.FontSize = 15; 
xticks(0.01:0.01:0.07);
xlabel('\alpha');
ylabel('BIC');
% plot(alpha_list,r1,'r-*');
% plot(alpha_list,r2,'b-^');
plot(alpha_list,r3,'k-x');
grid minor
% legend({'r = 1','r = 2','r = 3'});
legend({'r = 3'});
legend('Location','Best');
hold off
    
    
    