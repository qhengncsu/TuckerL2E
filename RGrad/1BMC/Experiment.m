function result = Experiment(pct,d,r,maxInner, sigmaVals, trials)
%% Performs synthetic experiments under probit (Gaussian) model
%
% Usage:  result = Experiment(pct,d,r,maxInner, sigmaVals, trials)
% Inputs:  pct - percentage of the entries in the matrix to observe (0-100)
%          d - dimension of synthetic d x d matrix to generate
%          r - rank of synthetic d x d matrix
%          maxInner - maximum number of iterations in spgSolver
%          sigmaVals - vector of values of sigma (the standard deviation of
%                      the noise in the Gaussian model)
%          trials - number of trials
% Outputs: result - has the following fields
%                   err   size(1) = # of sigma values
%                         size(2) = # of trials
%                         size(3) = nStat: # of error statistics per method
%
% Error statistics are defined/computed below in the code.
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


%% Setup for number of error statistics computed
nStat = 4;          % # of error statistics per method 
nErr  = 2* nStat;   % Two methods (1 = Combined, 2 = Nuclear norm)


%% Run separate experiments for each sigma and combine the results
% ------------------------------------------------------------------------
if length(sigmaVals) > 1
   n = length(sigmaVals);
   
   % Pre-allocate result arrays
   result = struct();
   result.err    = zeros(n,trials,nErr);
   result.normM  = zeros(n,trials);
   result.obsErr = zeros(n,trials);
   result.nIter  = zeros(n,trials,2);
   result.iSolve = cell(n,trials,2);
   
   for i=1:length(sigmaVals)
      data = Experiment(pct,d,r,maxInner,sigmaVals(i), trials);
      result.err(i,:,:)  = data.err;
      result.normM(i,:)  = data.normM;
      result.obsErr(i,:) = data.obsErr;
      result.nIter(i,:,:)  = data.nIter;
      result.iSolve(i,:,:) = data.iSolve;
   end

   return ;
end


%% Single sigma level from here
% ------------------------------------------------------------------------
d1   = d;
d2   = d;

%% Create file name for storing results
% ------------------------------------------------------------------------
p        = mfilename('fullpath');                         
prefix   = fileparts(p);                          % Base path (cur. dir.)
prefix   = [prefix,filesep,'cache'];              % Cache path
sigma    = power(10,sigmaVals);
str      = sprintf('%03d',round(sigmaVals*100));  % Sigma value
str      = strrep(str,'-','n');                   % Replace negative by n
filename = sprintf('%s/Exp_%03d_%d_%d_%d_%s.mat',prefix,pct,d,r,maxInner,str);


%% Load existing data when available, otherwise initialize new data
% ------------------------------------------------------------------------
if (exist(filename,'file'))
   data = load(filename);
   
   % Load random number generator
   s = data.s;
   
   err    = data.err;    % Various error metrics
   normM  = data.normM;  % Norm of M
   obsErr = data.obsErr; % Number of "observation errors"
   nIter  = data.nIter;  % Solver iterations
   iSolve = data.iSolve; % Solver info
   
   % Allocate extra memory for the arrays if needed
   if (trials > data.trial)
      err(trials,nErr) = 0;
      normM(trials,1)  = 0;
      obsErr(trials,1) = 0;
      nIter(trials,1)  = 0;
      iSolve{trials,1} = [];
   end

   % Start trial
   trialStart = data.trial + 1;
else
   % Initialize random number generator
   s = RandStream.create('mt19937ar','seed',5488);

   err    = zeros(trials,8);
   normM  = zeros(trials,1);
   obsErr = zeros(trials,1);
   nIter  = zeros(trials,2);  % Solver iterations -- combined & NN
   iSolve = cell(trials,2);   % Solver info       -- combined & NN

   trialStart = 1;
end


%% Run remaining trials
% ------------------------------------------------------------------------
for j = trialStart:trials
    
    % Generate a random low-rank matrix using uniform matrices
    M1 = rand(s,d1,r)-.5;
    M2 = rand(s,d2,r)-.5;
    M = M1*M2';    
    M = M/max(abs(M(:)));
    normM(j,1) = norm(M,'fro');
   

    % ------ Single noise level ------
    
    % Define observation model
    f      = @(x) gausscdf(x,0,sigma);
    fprime = @(x) gausspdf(x,0,sigma);

    % Obtain signs of noisy measurements
    Y = sign(f(M)-rand(s,d1,d2));
    y = Y(:);
    
    % Observe 'pct' percent of the entries in Y
    idx = find(rand(s,d1*d2,1) <= pct/100);
    m = length(idx);
    
    % Calculate number of "observation errors" as a fraction of m
    obsErr(j,1) = length(find(y(idx)~= sign(f(M(idx))-.5)))/m;

    % Set up and solve optimization problem
    options = struct();
    options.iterations = maxInner; % 10000
    options.stepMax    = 10000;
    options.stepMin    = 1e-4;
    options.optTol     = 1e-3;
    options.stepMax    = 1e9; % 1e9 1e8;
        
    funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);

    % Define alpha to be the correct maximum using an oracle
    alpha   = 1;
    radius  = alpha * sqrt(d1*d2*r);

    % Run both methods (1 = Combined nuclear and infinity norm, 2 = Nuclear norm only)
    for method=1:2
       if (method == 1)
          funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);
       else
          funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);
       end

       % Call solver
       [xhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
       xhat = reshape(xhat,d1,d2);

       [U,S,V] = svd(xhat);
       x1_debias = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';

       % Error statistics
       v = zeros(1,nStat);
       v(1) = norm(M-xhat,'fro');                    % Frobenius norm
       v(2) = norm(M-x1_debias,'fro');               % Frobenius norm (after debias)
       v(3) = sum(sign(xhat(:)) ~= sign(M(:)));      % Number of incorrect signs
               
       % Hellinger distance
       fM = f(M);
       fMhat = f(xhat);
       H = (sqrt(fM) - sqrt(fMhat)).^2 + (sqrt(1-fM) - sqrt(1-fMhat)).^2;
       v(4) = sum(H(:))/(d1*d2);    % Hellinger distance
           
       % Copy v vector into the error statistics
       err(j,(method-1)*nStat + (1:nStat)) = v;
    
       % Add objective value at x
       [fValM,gValM] = funObj(M(:));
       info.fValM  = fValM;
       info.gNormM = norm(gValM,2);
       [U,S,V] = svd(M);
       info.normInfM = max(abs(M(:)));
       info.normNucM = sum(abs(diag(S)));
    
       % Collect solver information
       nIter(j,method)  = info.iter;
       iSolve{j,method} = info;
       
    end % Methods
    
    % Save intermediate result
    trial = j;
    save(filename,'err','normM','obsErr','nIter','iSolve','s','trial');
end

%% Prepare result
% ------------------------------------------------------------------------
result = struct();
result.err    = err;
result.normM  = normM;
result.obsErr = obsErr;
result.nIter  = nIter;
result.iSolve = iSolve;
