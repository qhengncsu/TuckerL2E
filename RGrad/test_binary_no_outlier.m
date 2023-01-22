clc;clear;
addpath(genpath('functions'));

% fixed parameters
param.d = [100,100,100];
param.r = [2,2,2];

% % not fixed parameters
% beta_list = [0.1,0.5,1,2,5,10]; % stepsize
% mu0_list = [2,3,5,10,20];
% zeta_list = [1,5,10,20];
% sigma_list = [0.1,0.5,1,2,5];
% perturb_sigma_list = [1,2,5,10,15];

%%
param.mu0 = 5;

param.beta = 2.2;
param.sigma = 0.5;

param.zeta = 30;
param.perturb_sigma = 25;

param.init_hosvd = 0;

binary_no_outlier(param);