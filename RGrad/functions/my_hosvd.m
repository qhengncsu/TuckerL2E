function [t,nC,nU,nUT] = my_hosvd(C,U,UT,G,V)
% returns hosvd(T - P_{T}G, r), where T = tprod(C,U) and C is of low
% dimension, r = size(C)
% V is for calculating P_{T}G

r = size(C);
m = length(U);
L = zeros(2*r);
D = tprod(G,UT);

% need to modify this part for general m
L(1:r(1),1:r(2),1:r(3)) = C - D;
L(r(1)+1:2*r(1),1:r(2),1:r(3)) = -C;
L(1:r(1),r(2)+1:2*r(2),1:r(3)) = -C;
L(1:r(1),1:r(2),r(3)+1:2*r(3)) = -C;

Q = cell(m,1);
R = cell(m,1);

for i = 1:m
    X = [U{i} V{i}];
    [Q{i},R{i}] = qr(X,0);
end

tL = tprod(L,R);
tU = cell(m,1);
tUT = cell(m,1);

for i = 1:m
    X = ndim_unfold(tL,i);
    [uu,~,~] = svd(X,'econ');
    ur = uu(:,1:r(i));
    tU{i} = ur;
    tUT{i} = ur';
end

nU = cell(m,1);
nUT = cell(m,1);
for i = 1:m
    nU{i} = Q{i}*tU{i};
    nUT{i} = nU{i}';
end
nC = tprod(tL,tUT);
t = tprod(nC,nU);




