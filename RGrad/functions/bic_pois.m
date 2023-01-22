function v = bic_pois(A,S,T,r,I)

d = size(A);
dprod = prod(d);
v = (nnz(S) + sum(d.*r))*log(dprod) - (2*sum(A.*(T+S),'all') - 2*I*sum(exp(T+S),'all'))/I;