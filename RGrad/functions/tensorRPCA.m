function result = tensorRPCA(argin)
maxit = 100;

alpha = argin.alpha;
gamma = argin.gamma;
beta = argin.beta;
mu0 = argin.mu0;
d = argin.d;
r = argin.r; % rank of the tensor
sigma = argin.sigma;
stop_thres = argin.stop_thres;

%% generate A
rng(1);
Ttrue = randn(d);
[Ttrue,~,~,~,~,~] = hosvd(Ttrue,r); 
Ttrue = trimtensor(Ttrue,mu0); 

rng(alpha*100000+17);
Strue = binornd(ones(d),alpha) .* randn(d);

rng(sigma*100000+77);
Z = randn(d)*sigma;
A = Ttrue + Strue + Z;

%% obtain an initialization
relerrList = zeros(maxit,1);
% relerrSList = zeros(maxit,1);
relerrTinftyList = zeros(maxit,1);
% relerrSinftyList = zeros(maxit,1);

init_hosvd = argin.init_hosvd;
if init_hosvd == 0
    perturb_sigma = argin.perturb_sigma;
    rng(2);
    perturb = randn(d)*perturb_sigma;
    T = Ttrue + perturb; 
    [T,C,U,UT,~,~] = hosvd(T,r); 
else
    [T,C,U,UT,~,~] = hosvd(A,r); 
end
    
S = threshold(T,gamma,alpha,A); 

init_err = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
disp(['initialization rel err.: ', num2str(init_err)])
% disp(['initialization rel err. (S): ', num2str(norm(S(:) - Strue(:))/norm(Strue(:)))])

%% iteration
tic;
Tprev = T;
for i = 1:maxit
    G = T + S - A;
    [G_p,~,~] = mani_proj(G,C,U,UT);
    W = T - beta*G_p;
    [T,C,U,UT,~,~] = hosvd(W,r); 
    tildeT = trimtensor(T,mu0); % norm(tildeT(:))
    S = threshold(tildeT,gamma,alpha,A); 
    
    relerrList(i,1) = norm(tildeT(:) - Ttrue(:))/norm(Ttrue(:)); 
    relerrTinftyList(i,1) = norm(tildeT(:) - Ttrue(:), inf)/norm(Ttrue(:),inf);
%     relerrSList(i,1) = norm(S(:) - Strue(:))/norm(Strue(:));
%     relerrSinftyList(i,1) = norm(S(:) - Strue(:),inf)/norm(Strue(:),inf);
    
    disp(['iter:', num2str(i),', rel err.:',num2str(relerrList(i,1))])
    if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres
        break;
    end
    Tprev = tildeT;
end
runtime = toc

%% result
result.runtime = runtime;
result.param = argin;
result.relerrList = relerrList;
result.init_err = init_err;
result.actual_iter = i;
% result.outT = tildeT;
% result.outS = S;
result.relerrTinftyList = relerrTinftyList;
% result.relerrSList = relerrSList;
% result.relerrSinftyList = relerrSinftyList;



