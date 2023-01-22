function [T,C,U,UT] = HOOI(A,r)

maxit = 100;
sz = size(A);
m = length(sz);

U = cell(m,1);
UT = cell(m,1);

for i = 1:m
    U{i} = svdtrunc(ndim_unfold(A,i),r(i));
    UT{i} = U{i}';
end

for t = 1:maxit
    for i = 1:m
        UT{i} = eye(sz(i));
        X = tprod(A,UT);
        U{i} = svdtrunc(ndim_unfold(X,i),r(i));
        UT{i} = U{i}';
    end
end

C = tprod(A,UT);
T = tprod(C,U);