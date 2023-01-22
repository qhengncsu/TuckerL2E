clc;clear;
addpath(genpath(pwd));
%%
load 'data/it4tnsr_notime.mat';
A = log(1+it4tnsr_notime);

rng('default');
per = 0.1;
[A,testA,nnz_idx] = data_partition(A,per); % A for estimating T* and S*

%%
r = [3,3,3];
alpha = 0.003;
gamma = 1;

[T,S,time] = tRPCA_RGrad(A,r,gamma,alpha);

err = (T+S).*nnz_idx - testA;
norm(err(:))
