function movielens(mu)
%% Performs experiment on movielens dataset
%
% This file contains the code for reproducing the simulations in the paper 
% http://arxiv.org/abs/1209.3672 on the MovieLens dataset. For the precise 
% parameters used in the paper, the results are cached and are 
% automatically retrieved from the "cache" folder. If you modify these 
% parameters, the new results will be added to the "cache" folder.
%
% This code assumes that you have installed TFOCS (available at
% http://cvxr.com/tfocs/) and that the directory containing the function 
% "solver_sNuclearBPDN" has been added to your path.  It also assumes that
% the MovieLens 100k dataset have been loaded into a variable called
% "ratings" and saved in the file "ratings_100k.mat".  This can be created
% by downloading the ml-100k.zip file from http://www.grouplens.org/node/73
% and then extracting "u.data" and running the commands
% >> x = load('u.data');
% >> ratings = x(:,1:3);
% >> save ratings_100k ratings
%
% Usage:  movielens(mu)
% Inputs:  mu - regularization parameter for TFOCS
%
% Most recent change - 10/23/2013
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

s = RandStream.create('mt19937ar','seed',5489);

addpath TFOCS % Adjust if TFOCS is installed elsewhere
warning off

% Load 100k Movie Lens dataset
load('ratings_100k')
M = spconvert(ratings);
[d1 d2] = size(M);
numRatings = length(ratings(:,1));

% Generate observed data by comparing to average rating in dataset
averageRating = mean(ratings(:,3));
Y = sign(M-averageRating);
y = Y(:);

% Observe 95 percent of the sampled entries in Y
subIdx   = find(rand(s,numRatings,1) < .95); 
trainIdx = sub2ind(size(M),ratings(subIdx,1),ratings(subIdx,2)); % Convert to vectorized indices

% Define 1-bit and traditional (Frobenius norm) objective functions
f       = @(x) (1 ./ (1 + exp(-x)));
fprime  = @(x) (exp(x) ./ (1 + exp(x)).^2);     
funObj1 = @(x) logObjectiveGeneral(x,y,trainIdx,f,fprime);       

% Set options for SPG solver
options = struct();
options.iterations = 1000; 
options.stepMax    = 10000;
options.stepMin    = 1e-4;
options.optTol     = 1e-3;
options.stepMax    = 1e9; % 1e9 1e8;
options.verbosity  = 1;

% Iterate over values for nuclear norm constraint parameter
numScales = 20;
scaleVals = logspace(-.5,1,numScales);
froVals = (norm(M(trainIdx))/2)*logspace(-5,-.5,numScales);

% Pre-allocate results
spgOutput   = cell(1,numScales);
tfocsOutput = cell(1,numScales);
Pcorrect1   = zeros(1,numScales);
Pcorrect2   = zeros(1,numScales);
Pcorrect3   = zeros(1,numScales);
Pcorrect4   = zeros(1,numScales);
PeMat1      = zeros(numScales,5);
PeMat2      = zeros(numScales,5);
PeMat3      = zeros(numScales,5);
PeMat4      = zeros(numScales,5);

filename = sprintf('cache/movielensResults_%d.mat',round(100*mu));

if (exist(filename)) 
   error('Results already stored in cache.');
end

for jj = 1:numScales,
   
    fprintf('==========================================================\n');
    fprintf('%d\n',jj);
    fprintf('==========================================================\n');
    

    % Set radius and define projection function
    radius  = scaleVals(jj)*sqrt(d1*d2);
    funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);

    % 1-bit matrix completion
    [Mhat1,info1] = spgSolver(funObj1, funProj, zeros(d1*d2,1), options);
    % Traditional matrix completion
    [Mhat2,info2] = solver_sNuclearBPDN( {d1 d2 trainIdx}, M(trainIdx), froVals(jj), mu);
    info2 = rmfield(info2,'dual');
    
    spgOutput{jj} = info1;
    tfocsOutput{jj} = info2;
    
    
    Yhat1 = sign(f(Mhat1)-.5);
    Yhat2 = sign(Mhat2-averageRating);
    
    v3 = Mhat2(trainIdx);
    Yhat3 = sign(Mhat2-mean(v3(:)));
    Yhat4 = sign(Mhat2-mean(Mhat2(:)));
    
    subIdx2 = setdiff([1:numRatings],subIdx);
    testIdx = sub2ind(size(M),ratings(subIdx2,1),ratings(subIdx2,2));
    
    err = y(testIdx)-Yhat1(testIdx);
    Pcorrect1(jj) = 1-length(find(err~=0))/length(err);
    
    err = y(testIdx)-Yhat2(testIdx);
    Pcorrect2(jj) = 1-length(find(err~=0))/length(err);

    err = y(testIdx)-Yhat3(testIdx);
    Pcorrect3(jj) = 1-length(find(err~=0))/length(err);

    err = y(testIdx)-Yhat4(testIdx);
    Pcorrect4(jj) = 1-length(find(err~=0))/length(err);

    for mm=1:5,
        subIdx3 = find(M==mm);
        testIdxSub = intersect(subIdx3,testIdx);
        err = y(testIdxSub)-Yhat1(testIdxSub);
        PeMat1(jj,mm) = 1-length(find(err~=0))/length(err);
        
        err = y(testIdxSub)-Yhat2(testIdxSub);
        PeMat2(jj,mm) = 1-length(find(err~=0))/length(err);

        err = y(testIdxSub)-Yhat3(testIdxSub);
        PeMat3(jj,mm) = 1-length(find(err~=0))/length(err);

        err = y(testIdxSub)-Yhat4(testIdxSub);
        PeMat4(jj,mm) = 1-length(find(err~=0))/length(err);
    end
    
    disp(['Results for iteration ' num2str(jj) ' : ' num2str(Pcorrect1(jj)) ' -- ' num2str(Pcorrect2(jj))])

    save(sprintf('cache/movielensResults_%d.mat',round(100*mu)),'spgOutput','tfocsOutput',...
         'Pcorrect1','Pcorrect2','Pcorrect3','Pcorrect4',...
         'PeMat1','PeMat2','PeMat3','PeMat4');
end

