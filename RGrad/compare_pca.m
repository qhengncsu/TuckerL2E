clc;clear;
addpath(genpath('functions'));
%%
d = [300,300,300];
r = [2,2,2];

beta = 0.5; % stepsize
mu0 = 5; % incoherence parameter
% sigma = 0.01; % noise level
sigma_set = 0.01:0.01:0.05;
tol = 0.001;
repeat = 10;

alg = 2; % 1 -- RGrad; 2 -- PGD

%% 

result = zeros(length(sigma_set),repeat);
total_time_list = zeros(length(sigma_set),1);
avg_time_list = zeros(length(sigma_set),1);

%% generate the observed tensor A
rng('default');
for k = 1:length(sigma_set)
    total_time = 0;
    total_iter = 0;
    sigma = sigma_set(k);
    for tt = 1:repeat
        Ttrue = 10 * randn(d);
        [Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
        Ttrue = trimtensor(Ttrue,mu0); 
        
        Z = randn(d)*sigma;
        A = Ttrue + Z;
        fprintf('=====================================================\n');
        fprintf('Current noise level = %f, repeating %d times\n',sigma,tt);
        %% initialization 
        [T,C,U,UT,~,~] = hosvd(A,r); % use HOSVD as initialization
%         [T,C,U,UT] = initialization(A,r);
        init_err = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        disp(['initialization rel err.: ', num2str(init_err)])
        %% iteration
        maxit = 100;

        if alg == 1 %RGrad
            tic;
            for i = 1:maxit
                Tprev = T;
                G = T - A;
                [G,~,V] = mani_proj(G,C,U,UT);
                [T,C,U,UT] = my_hosvd(C,U,UT,beta*G,V);
                disp(['iter:', num2str(i),', rel err.:',num2str(norm(T(:) - Ttrue(:))/norm(Ttrue(:)))])
                if norm(T(:) - Tprev(:))/norm(T(:)) < tol
                    break;
                end
            end
            time = toc;
            total_time = total_time + time;
            total_iter = total_iter + i;
            result(k,tt) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
        elseif alg == 2 % PGD
            tic;
            for i = 1:maxit
                Tprev = T;
                G = T - A;
                W = T - beta*G;
                [T,C,U,UT,~,~] = hosvd(W,r); 

                disp(['iter:', num2str(i),', rel err.:',num2str(norm(T(:) - Ttrue(:))/norm(Ttrue(:)))])
                if norm(T(:) - Tprev(:))/norm(T(:)) < tol
                    break;
                end
            end
            time = toc;
            total_time = total_time + time;
            total_iter = total_iter + i;
            result(k,tt) = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
            
        end
    end
    total_time_list(k,1) = total_time;
    avg_time = total_time/total_iter;
    avg_time_list(k,1) = avg_time;
end

%%
myResult.result = result;
myResult.total_time_list = total_time_list;
myResult.avg_time_list = avg_time_list;
myResult.alg = alg;
save('RGrad_result_10_1.mat','myResult')
%%
% fprintf('||Z||_F/||T*||_F = %f\n',norm(Z(:))/norm(Ttrue(:)));
% fprintf('Total runtime = %f; Average runtime = %f\n',time,time/i);
% fprintf('Relative error = %f\n',norm(T(:) - Ttrue(:))/norm(Ttrue(:)));
    
    

