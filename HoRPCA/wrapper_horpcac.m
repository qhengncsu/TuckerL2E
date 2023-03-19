function results = wrapper_HoRPCAC(X,W,R)
data.X = X;
if isa(W,'tensor')
  data.linInd = find(W(:));
  data.b = X(data.linInd);
else
  data.b = X(:);
  data.linInd = reshape(1:prod(size(X)),prod(size(X)),1);
end
params.X0 = tenzeros(size(data.X));
N = ndims(data.X);
params.V0= cell(1,N);
for i = 1:N
    params.V0{i} = tenzeros(size(data.X));
end
params.E0 = tenzeros(size(data.X));
params.mu0 = 1/(N+1);
%T = double(data.T);
if ~exist( 'mu1fac', 'var' ) || isempty(mu1fac); mu1fac = 5; end
params.mu1fac = mu1fac;
params.mu1 = mu1fac*std(data.b);  %mu1fac*std(T(:));
params.mu2 = params.mu1;
params.mu_min = 1e-4;
params.mu_max = 1e2;
params.max_iter = 1000;
params.opt_tol = 1e-3;  %1e-5 for syn data analysis plots, 1e-3 for others
params.eta = 1/(N+1);
params.IsTC = exist( 'IsTC', 'var' ) && ~isempty(IsTC) && IsTC;
if isscalar(R)
  R = ones(1,N)*R;
end
params.k = R;

r = 1/sqrt( max(size(data.X)) ); %0.015; %0.09:-0.005:0.02;
if ~exist( 'lambdaS', 'var' ) || isempty(lambdaS)
    lambdaS = 1;
end
params.lambdaS = lambdaS;     %1e1;
if ~exist( 'rRatio', 'var' ) || isempty(rRatio)
    rRatio = 1/4;
end
params.lambda = params.lambdaS*r*rRatio;       %params.lambdaS*r/4;
% params.lambda = 0.06;
params.rRatio = rRatio;
params.verbose = exist( 'verbose', 'var' ) && ~isempty(verbose) && verbose;
params.use_cont = true;
%if ~exist( 'mode', 'var' ) || isempty(mode)
    %mode = N;
%end
%%%%%%%%%% for PROPACK %%%%%%%%%%%%
% declare global var 'sv'
global sv;
global tmode;
global use_propack;
global curr_mu;
sv =  ceil(min(size(data.X)) * 0.1) * ones( 1, N );
use_propack = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = tensor_rpca_tc_adal_ncx(data,params);
end