function trimT = trimtensor(T,mu0)

sz = size(T);
m = length(sz);
r = zeros(m,1);
for i = 1:m
    r(i) = rank(ndim_unfold(T,i));
end
[~,C,U,~,~,~] = hosvd(T,r);
trimU = cell(m,1);

for i = 1 : length(U)
    trimU{i} = trim(U{i},mu0);
end

trimT = tprod(C,trimU);