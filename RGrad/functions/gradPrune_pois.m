function S = gradPrune_pois(T,alpha,kpr,A,I)

G = -A/I + exp(T);
J = findAInd(G,alpha);
S = log(A/I) - T;
S(S>kpr) = kpr;
S(S<-kpr) = -kpr;
S = S .* J;