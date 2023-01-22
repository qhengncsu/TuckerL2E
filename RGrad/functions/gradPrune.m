function S = gradPrune(T,alpha,kpr,A)

G = -A./(1+exp(T)) + (1-A)./(1+exp(-T));
J = findAInd(G,alpha);
S = (2*A - 1)*kpr - T;
S = S .* J;

