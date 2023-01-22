function T = generate_low_rank(d,r)

C = randn(r)-0.5;
m = length(d);
U = cell(m,1);
for i = 1:m
    U{i} = rand(d(i),r(i))-0.5;
    [q,~] = qr(U{i},0);
    U{i} = q;
end

T = tprod(C,U);
T = T/(max(abs(T(:))));