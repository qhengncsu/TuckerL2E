function [Y,tau,out]=tucker_l2e_opt(Z,R,varargin)
%% Error checking
if ~isa(Z,'tensor') && ~isa(Z,'sptensor')
  error('Z must be a tensor or a sptensor');
end
 
if (nargin < 2)
    error('Error: invalid input arguments');
end

sz = size(Z);
N = length(sz);

params = inputParser;
params.addParameter('W',"all_observed", @(x) isstring(x)||isa(x,'tensor'));
params.addParameter('init', 'tucker_als', @(x) ismember(x,{'spectral','tucker_als','HOSVD'}));
params.addParameter('opt_options', '', @isstruct);
params.addParameter('scale', 0.1, @(x) isnumeric(x));
params.addParameter('taumax', 50, @(x) isnumeric(x));
params.addParameter('maxiters', 5000, @(x) isnumeric(x));

params.parse(varargin{:});

W = params.Results.W;
init = params.Results.init;
options = params.Results.opt_options;
s = params.Results.scale;
taumax = params.Results.taumax;
maxiters = params.Results.maxiters;
if isempty(options) 
  options.factr = 1e7;
  options.m = 3;
  options.maxIts = maxiters;
  options.maxTotalIts = 500000;
  options.printEvery = 100;
end

if isscalar(R)
  R = ones(1,N)*R;
end

coreparams = prod(R);
facparams = sum(sz.*R);
Nparams = coreparams+facparams+1;
l = -inf(Nparams,1);
u = inf(Nparams,1);

if isa(W,'tensor')
  Z = Z.*W;
end

if isa(W,'tensor')
  total_entries = collapse(W);
else
  total_entries = prod(sz);
end

Z_avg = collapse(Z)/total_entries;

data = Z(:);
nonzero_entries = nonzeros(data);
maddata = mean(abs(nonzero_entries-Z_avg),'all');
scale = maddata/s

Z = Z/scale;
Z_avg = Z_avg/scale;
if init=="spectral"
  p = total_entries/prod(sz);
  A = cell(N,1);
  for i=1:N
      T = double(tenmat(Z,i))/p;
      gram = T*(T.');
      gram = gram - diag(diag(gram));
      [U,~,~] = svds(gram,R(i));
      A{i} = U;
  end
  G = ttm(Z./p,A,'t');
elseif init=="tucker_als"
  if isa(W,'tensor')
    Z = Z + (1-W).*Z_avg;
  end
  Y = tucker_als(Z,R);
  G = Y.core;
  A = Y.U;
elseif init=="HOSVD"
  if isa(W,'tensor')
    Z = Z + (1-W).*Z_avg;
  end
  Y = hosvd_tt(Z,1,'ranks',R);
  G = Y.core;
  A = Y.U;
end

tau = 0.01;

function x = fac_to_vec(G,A,logtau)
x = zeros(Nparams,1);
x(1:coreparams) = G(:);
for n = 1:N
  idx1 = coreparams+sum(sz(1:n-1).*R(1:n-1))+1;
  idx2 = coreparams+sum(sz(1:n).*R(1:n));
  x(idx1:idx2) = reshape(A{n},sz(n)*R(n),1);
end
x(end) = logtau;
end

function [G,A,logtau] = vec_to_fac(x)
G = tensor(x(1:coreparams),R);
A = cell(N,1);
for n = 1:N
  idx1 = coreparams+sum(sz(1:n-1).*R(1:n-1))+1;
  idx2 = coreparams+sum(sz(1:n).*R(1:n));
  A{n} = reshape(x(idx1:idx2),sz(n),R(n));
end
logtau = x(end);
end

function [f,gG,gA,glogtau] = tucker_fg_GA(Z,G,A,W,logtau)
tau = exp(logtau);
L = ttm(G,A);
if isa(W,'tensor')
  Res = W.*(Z-L);
  E = tenfun(@(x)(exp(-x.^2.*tau^2/2)),Res).*W;
else
  Res = Z - L;
  E = tenfun(@(x)(exp(-x.^2.*tau^2/2)),Res);
end
f = total_entries/(2*sqrt(pi))*tau - tau*sqrt(2/pi)*collapse(E);

if nargout>1
  T = tenfun(@(x)(tau^2.*(x.^2)-1),Res);
  core = -tau^3*sqrt(2/pi)*E.*Res;
  ET = E.*T;

  gA = cell(N,1);
  for n=1:N
    gA{n} = double(tenmat(core,n))*transpose(double(tenmat(ttm(G,A,-n),n)));
  end
  gG = ttm(core,A,'t');

  glogtau = (total_entries/(2*sqrt(pi))+sqrt(2/pi)*collapse(ET))*tau;
end
end

function [f,g] = tucker_fun(x,Z,W)
[G,A,logtau] = vec_to_fac(x);
if nargout == 1
  f = tucker_fg_GA(Z,G,A,W,logtau);
else
  [f,gG,gA,glogtau] = tucker_fg_GA(Z,G,A,W,logtau);
  g = fac_to_vec(gG,gA,glogtau);
end
end

opts = options;
opts.x0 = fac_to_vec(G,A,log(tau));
u(end) = log(taumax);
[xx,ff,out] = lbfgsb(@(x)tucker_fun(x,Z,W), l, u, opts);

[G,A,logtau] = vec_to_fac(xx);
tau = exp(logtau);
G = G*scale;
Y = ttensor(G,A);
end