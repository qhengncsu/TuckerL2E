function [T,C,U,UT] = initialization(A,r)

tau = 1;
wtA = trunc(A,tau);
[T,C,U,UT] = HOOI(wtA,r);