%% syntheticSims.m
% 
% This file contains the code for reproducing the synthetic simulations in
% the paper http://arxiv.org/abs/1209.3672. For the precise parameters
% used in the paper, the results are cached and are automatically retrieved
% from the "cache" folder. If you modify these parameters, the new results
% will be added to the "cache" folder.
%
% Most recent change - 10/3/2013
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

%% Synthetic simulations for fixed number of observations, varying noise level
sigmaVals = [-3:.25:1]; % Various noise levels
n = 500; % Dimension of matrix
r = 1; % Rank

% Synethtic simulations observing 15% of entries and doing 15 trials
results = Experiment(15,n,r,10001,sigmaVals,15); 

%% Synthetic simulations for fixed noise level, varying number of observations
sigma = sigmaVals(10); % Fixed noise level
n = 200; % Dimension of matrix
mVals = [5:5:100]; % Various numbers of observations

mNum = length(mVals);    
results1 = cell(mNum,1);
results2 = results1;
results3 = results1;

r = 3; % Rank 3
for jj=1:mNum,
    results1{jj} = Experiment(mVals(jj),n,r,10001,sigma,15);
end

r = 5; % Rank 5
for jj=1:mNum,
    results2{jj} = Experiment(mVals(jj),n,r,10001,sigma,15);
end

r = 10; % Rank 10
for jj=1:mNum,
    results3{jj} = Experiment(mVals(jj),n,r,10001,sigma,15);
end


%% Plot Frobenius error as a function of sigma
figure()
plot(sigmaVals,mean((results.err(:,:,2)./results.normM)'),'b','Linewidth',2,'Marker','x')
hold on
plot(sigmaVals,mean((results.err(:,:,6)./results.normM)'),'r:','Linewidth',2,'Marker','o')

legend('Infinity-norm constraint','No infinity-norm constraint')


%% Plot Frobenius error as a function of m
for j=1:20,
    error1(j) = mean(results1{j}.err(:,6)./results1{j}.normM);
    error2(j) = mean(results2{j}.err(:,6)./results2{j}.normM);
    error3(j) = mean(results3{j}.err(:,6)./results3{j}.normM);
end

figure()
plot(mVals(3:end)/100,error1(3:end),'b','Linewidth',2,'Marker','x')
hold on
plot(mVals(3:end)/100,error2(3:end),'r:','Linewidth',2,'Marker','o')
plot(mVals(3:end)/100,error3(3:end),'g-.','Linewidth',2,'Marker','square')
axis([0 1 0 2])

legend('r=3','r=5','r=10')


%% Plot Hellinger distance as a function of m

for j=1:20,
    error1(j) = mean(results1{j}.err(:,8));
    error2(j) = mean(results2{j}.err(:,8));
    error3(j) = mean(results3{j}.err(:,8));
end

figure()
plot(mVals(3:end)/100,error1(3:end),'b','Linewidth',2,'Marker','x')
hold on
plot(mVals(3:end)/100,error2(3:end),'r:','Linewidth',2,'Marker','o')
plot(mVals(3:end)/100,error3(3:end),'g-.','Linewidth',2,'Marker','square')
axis([0 1 0 .5])

legend('r=3','r=5','r=10')

