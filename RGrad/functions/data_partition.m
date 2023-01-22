function [trainA,testA,nnz_idx] = data_partition(A,per)

sz = size(A);
nnz_idx = A~=0;
T = rand(sz);
omega = T<per;
omega = omega.*nnz_idx;

testA = A.*omega;
trainA = A - testA;

nnz_idx = omega;

