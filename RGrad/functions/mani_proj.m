function [ A_p, G, V ] = mani_proj(A, S, U, UT)
% Project gradient to tangent space 
% input: A - the tensor to be projected
%        S,U - T = tprod(S,U)
%        UT - transpose of U
% output: A_p = P_T(A)
%         G,V - intermediate variable
              

d = length(size(A));
G = tprod(A, UT);
A_p = tprod(G, U);
V = cell(1, d);
U_c = U;
UT_c = UT;
for i = 1:d
    UT_c{i} = eye(size(A, i));
    S_unfoldi = ndim_unfold(S, i);
    
    t1 = ndim_unfold(tprod(A, UT_c), i)*S_unfoldi'*inv(S_unfoldi*S_unfoldi');
    t2 = U{i}'*t1;
    V{i} = t1 - U{i}*t2;
    
    U_c{i} = V{i};
    A_p = A_p + tprod(S, U_c);
    U_c = U;
    UT_c = UT;
end