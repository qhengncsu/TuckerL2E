function [T,C,U,UT] = Trim2(W,r,zeta)

W(W>zeta) = zeta;
W(W<-zeta) = -zeta;
[T,C,U,UT,~,~] = hosvd(W,r);