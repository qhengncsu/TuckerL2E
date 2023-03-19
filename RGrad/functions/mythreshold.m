function S = mythreshold(T,gamma,alpha,A)
% for tensor RPCA with loss function L(X) = 0.5||X-A||_F^2

G = T - A;
J = findAInd_new(G,gamma*alpha);
S = A - T;
S = S .* J;



