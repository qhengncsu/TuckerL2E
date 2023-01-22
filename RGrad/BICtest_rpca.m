% testing BIC criterion
clc;clear;
addpath(genpath(pwd));
%% set parameter
d = [100,100,100];
dprod = prod(d);
r0 = [3,3,3];
alpha0 = 0.1;
mu0 = 5;
sigma = 0.005;
gamma = 1;
stop_thres = 1e-3;
beta = 0.5;
maxit = 100;

%% generate A
rng('default');
Ttrue = randn(d);
[Ttrue,~,~,~,~,~] = hosvd(Ttrue,r0); 
Ttrue = trimtensor(Ttrue,mu0); 
Ttrue = 0.1*Ttrue/max(abs(Ttrue(:)));

Strue = binornd(ones(d),alpha0) .* randn(d);
Strue = 4 * Strue/max(abs(Strue(:)));

Z = randn(d)*sigma;
A = Ttrue + Strue + Z;

%%
r_list = [4,5]; r_len = length(r_list);
alpha_list = 0.02:0.02:0.2; alpha_len = length(alpha_list);
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
        fprintf('(%4f,%d)         %.4e        %.4e                 %.4e           %.4e\n',alpha,r_list(i),...
                       result(i,j),norm(T(:)-Ttrue(:))/norm(Ttrue(:)),(sum(abs(S(:))) + sum(d.*r)),dprod*log(norm(A(:)-T(:)-S(:))^2))
    end
end

%% visualization
r1 = result(1,:);
r2 = result(2,:);
r3 = result(3,:);
r4 = result(4,:);
r5 = result(5,:);

alpha_list = 0.02:0.02:0.2; alpha_list=[alpha_list(1:2) 0.05 alpha_list(3:end)];

hold on
ax = gca;
ax.FontSize = 15; 
xticks(0.02:0.02:0.2);
xlabel('\alpha');
ylabel('BIC');
plot(alpha_list,r1,'r-*');
plot(alpha_list,r2,'g-^');
plot(alpha_list,r3,'b-x');
plot(alpha_list,r4,'m-+');
plot(alpha_list,r5,'k-.');


grid minor
legend({'r = 1','r = 2','r = 3','r = 4','r = 5'});
legend('Location','Best');
hold off
    
    
    