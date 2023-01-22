function v = bic_binary(A,S,T,r,sigma)

d = size(A);
dprod = prod(d);
p = @(x) 1./(1+exp(-x/sigma));
v = (nnz(S) + sum(d.*r))*log(dprod) - 2*(sum(A.*log(p(T+S)),'all') + sum((1-A).*log(1-p(T+S)),'all')); %bic
% v = (sum(abs(S(:))) + sum(d.*r))*log(dprod) -2*(sum(A.*log(p(T+S)),'all') + sum((1-A).*log(1-p(T+S)),'all')); %bic+l1




% v = (nnz(S) + sum(d.*r)) -2* (sum(A.*log(p(T+S)),'all') + sum((1-A).*log(1-p(T+S)),'all')); %aic
% v = (sum(abs(S(:))) + sum(d.*r)) -2* (sum(A.*log(p(T+S)),'all') +sum((1-A).*log(1-p(T+S)),'all')); %aic+l1


