function [T,S,time] = tRPCA_RGrad(A,r,gamma,alpha,Ttrue,opts)
% Input: 
%     A     -- the observed tensor A = Ttrue + Strue + Z, where Ttrue is of low Tucker rank,
%               S is the sparse outlier, Z is additive noise
%     r     -- Tucker rank of Ttrue
%     alpha -- the sparsity
%     Ttrue -- low Tucker rank tensor to be recovered, only used for synthetic case for
%               calculating relative error
%     opts  -- Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.maxit      -   maximum number of iterations
%           opts.mu0        -   incoherence parameter
%           opts.verbose    -   display rel. error or not
%           opts.beta       -   step size
% Output:
%       T       -    Low rank part
%       S       -    Sparse part
%       time    -    runtime


tol = 1e-3;
maxit = 100;
beta = 0.3;
mu0 = 5;
verbose = 0;
exist_Ttrue = 1;

if nargin < 5
    exist_Ttrue = 0;
end

if ~exist('opts', 'var')
    opts = [];
end   

if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'maxit');       maxit = opts.maxit;          end
if isfield(opts, 'beta');        beta = opts.beta;            end
if isfield(opts, 'mu0');         mu0 = opts.mu0;              end
if isfield(opts, 'verbose');     verbose = opts.verbose;      end


%% initialization 
[T,C,U,UT,~,~] = hosvd(A,r); % use HOSVD as initialization
% [T,C,U,UT] = initialization(A,r);
S = mythreshold(T,gamma,alpha,A); 
if exist_Ttrue
    init_err = norm(T(:) - Ttrue(:))/norm(Ttrue(:));
%     fprintf('Initialization rel. err.: %f\n', init_err);
end
tic;
%% iteration
for i = 1:maxit
    Tprev = T;
    G = T + S - A;
    [G,~,V] = mani_proj(G,C,U,UT);
    [T,C,U,UT] = my_hosvd(C,U,UT,beta*G,V);        
    T = trimtensor(T,mu0);
    S = mythreshold(T,gamma,alpha,A); 
    if verbose && exist_Ttrue && (mod(i,10)==0 || i == 1)
        disp(['iter:', num2str(i),', rel err.:',num2str(norm(T(:) - Ttrue(:))/norm(Ttrue(:)))])
    end
    if norm(T(:) - Tprev(:))/norm(T(:)) < tol
        break;
    end
end
time = toc;


