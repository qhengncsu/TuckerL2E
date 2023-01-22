function [TT S U UT sv tol] = hosvd(T, r, dim)
%HOSVD High Order SVD of a multidimensional array
%	[S U sv tol] = HOSVD(T)
%	[S U sv tol] = HOSVD(T, dim)
%	[S U sv tol] = HOSVD(T, dim, tol)
%	
%	T    - multidimensional array
%	dim  - bools: hosvd is applied in these dimensions (default: ones())
%	tol  - tolerance for each dimensions (default: eps * ones())
%	
%	S    - decomposed core tensor so that T==tprod(S, U)
%	U    - matrices for each dimension (U{n}==[] if dim(n) was 0)
%	sv   - n-mode singular values (or [] if dim(n) was 0)
%	tol  - larges dropped singular value (hosvd truncates small sv)
%
%	
%	See also TPROD, SVDTRUNC

M = size(T);
P = length(M);
if nargin < 2
	r = M;
end
if nargin < 3
	dim = ones(1,P);
end
% if numel(tol) == 1
% 	tol = tol*ones(1,P);
% end
U = cell(1,P);
UT = cell(1,P);
sv = cell(1,P);
for i = 1:P
	if dim(i)
		A = ndim_unfold(T, i);
		% SVD based reduction of the current dimension (i)
		[Ui svi toli] = svdtrunc(A, r(i));
		U{i} = Ui;
		UT{i} = Ui';
		sv{i} = svi;
		tol(i) = toli;
	else
		U{i} = [];
		UT{i} = [];
		sv{i} = [];
		tol(i) = 0;
	end
end
S = tprod(T, UT);
% size(S)
TT = tprod(S, U);