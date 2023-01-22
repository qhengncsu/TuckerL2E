function v = bic(A,S,T,r)

d = size(A);
dprod = prod(d);
v = (nnz(S) + sum(d.*r))*log(dprod) + dprod*log(norm(A(:)-T(:)-S(:))^2); %bic
% v = (nnz(S) + sum(d.*r)) + dprod*log(norm(A(:)-T(:)-S(:))^2); %aic
% v = (sum(abs(S(:))) + sum(d.*r))*log(dprod) + dprod*log(norm(A(:)-T(:)-S(:))^2); %bic/l1
% v = sum(abs(S(:))) + sum(d.*r) + dprod*log(norm(A(:)-T(:)-S(:))^2); %aic/l1