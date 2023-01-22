function [T,C,U,UT] = init_binary(A,r,sigma)
sz = size(A);
Amat = ndim_unfold(A,1);
[d1,d2] = size(Amat);
rr = min(r(1),r(2)*r(3));
y = 2*Amat(:)-1;
idx = 1:length(y);

addpath('1BMC');
% Logistic model
f       = @(x) (1 ./ (1 + exp(-x/sigma)));
fprime  = @(x) (1/sigma) * (exp(x/sigma) ./ (1 + exp(x/sigma)).^2);  

%% Set up optimization problem
options = struct();
options.iterations = 10000; 
options.stepMax    = 10000;
options.stepMin    = 1e-4;
options.optTol     = 1e-3;
options.stepMax    = 1e9;
        
funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);

%% Define alpha to be the correct maximum using an oracle
alpha   = 1;
radius  = alpha * sqrt(d1*d2*rr);

%% Define constraints
% Use nuclear-norm constraint only
funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);

% Use nuclear-norm plus infinity-norm constraints
%funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);

%% Recover estimate Mhat of M
[Mhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
Mhat = reshape(Mhat,d1,d2);

[U,S,V] = svd(Mhat);
Mhat_debias = U(:,1:rr)*S(1:rr,1:rr)*V(:,1:rr)'; % Project onto actual rank if known

T = ndim_fold(Mhat,1,sz);
[T,C,U,UT] = Trim2(T,r,1);


