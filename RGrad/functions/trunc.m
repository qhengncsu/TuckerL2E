function T = trunc(A,tau)

T = A;
T(T > tau) = tau;
T(T < -tau) = -tau;