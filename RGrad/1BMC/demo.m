%% demo.m
% 
% This file contains the code for a simple demo of the one bit matrix 
% completion software package that supplements the paper 
% http://arxiv.org/abs/1209.3672. 
%
% Most recent change - 5/16/2014
%
% Copyright 2013, M. Davenport, Y. Plan, E. van den Berg, M. Wootters
%
% This file is part of 1BMC Toolbox version 1.2.
%
%    The 1BMC Toolbox is free software: you can redistribute it and/or 
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    The 1BMC Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the 1BMC Toolbox. If not, see <http://www.gnu.org/licenses/>.


%% Initialize
clear all;
close all;

d1 = 1e4; % Number of rows in matrix
d2 = 1e2; % Number of columns in matrix
r = 2; % Rank
pct = 100; % Percent of entries to observe
sigma = 3; % Noise level

strm = RandStream('mt19937ar','Seed',5);

%% Generate a random low-rank matrix using uniform matrices
M1 = rand(strm,d1,r)-.5;
M2 = rand(strm,d2,r)-.5;
M = M1*M2';    
M = M/max(abs(M(:)));

%% Define observation model (probit/Gaussian noise)
f      = @(x) gausscdf(x,0,sigma);
fprime = @(x) gausspdf(x,0,sigma);

% Logistic model
% f       = @(x) (1 ./ (1 + exp(-x)));
% fprime  = @(x) (exp(x) ./ (1 + exp(x)).^2);  

%% Obtain signs of noisy measurements
Y = sign(f(M)-rand(strm,d1,d2));
y = Y(:);
    
%% Observe 'pct' percent of the entries in Y
idx = find(rand(strm,d1*d2,1) <= pct/100);
m = length(idx);

%% Set up optimization problem
options = struct();
options.iterations = 10000; 
options.stepMax    = 10000;
options.stepMin    = 1e-4;
options.optTol     = 1e-3;
options.stepMax    = 1e9;
        
funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);

%% Define alpha to be the correct maximum using an oracle
alpha   = 1;
radius  = alpha * sqrt(d1*d2*r);

%% Define constraints
% Use nuclear-norm constraint only
funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);

% Use nuclear-norm plus infinity-norm constraints
%funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);

%% Recover estimate Mhat of M
[Mhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
Mhat = reshape(Mhat,d1,d2);

[U,S,V] = svd(Mhat);
Mhat_debias = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Project onto actual rank if known

%% Compute relative error
norm(Mhat_debias-M,'fro')/norm(M,'fro')
