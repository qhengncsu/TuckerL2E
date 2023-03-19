addpath('./HoRPCA/inexact_alm_rpca')
addpath('./HoRPCA/data')
addpath('./HoRPCA/lightspeed')
addpath('./HoRPCA/PROPACK')
addpath('./HoRPCA/code/rpca')
addpath('./HoRPCA/code/tc')
addpath('./HoRPCA/code/utils')
addpath('./HoRPCA')
addpath('./TuckerL2E')
addpath('./TuckerL2E/L-BFGS-B-C/Matlab')
addpath('./TuckerL2E/tensor_toolbox-v3.2.1')
addpath('./RGrad')
addpath('./RGrad/functions')
sz = [50 50 50];
reldiffs1_cpopt = zeros(9,1);
reldiffs1_tuckerals = zeros(9,1);
reldiffs1_horpcac = zeros(9,1);
reldiffs1_l2e = zeros(9,1);
reldiffs1_l2e_50 = zeros(9,1);
for j=1:1
  rng(j,'philox')
  info = create_problem('Type','CP','Size',sz,'Num_Factors',15,'Factor_Generator','randn','Noise',0);
  X = info.Data;
  Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
  truth = tensor(info.Soln);
  stdtruth = std(truth(:));
  S = tensor(zeros(sz));
  S(Oomega) = 10*stdtruth*(2*rand(length(Oomega),1)-1);
  X = X+S; 
  for i=1:9
    r = i*5;
    Xhat_cpopt = cp_opt(X,r,'init','nvecs');
    reldiffs1_cpopt(i,j) = norm(tensor(Xhat_cpopt)-truth)/norm(truth);
    Xhat_tuckerals = tucker_als(X,r);
    reldiffs1_tuckerals(i,j) = norm(tensor(Xhat_tuckerals)-truth)/norm(truth);
    result_horpcac = wrapper_horpcac(X,'all_observed',r);
    Xhat_horpcac = result_horpcac.X;
    reldiffs1_horpcac(i,j) = norm(tensor(Xhat_horpcac)-truth)/norm(truth);
    Xhat_l2e = tucker_l2e_opt(X,r,'taumax',50);
    reldiffs1_l2e(i,j) = norm(tensor(Xhat_l2e)-truth)/norm(truth);
  end
end    

cv = cvpartition(prod(sz),'KFold',10);
MADs_horpcac1 = zeros(9,1);
MADs_l2e1 = zeros(9,1);
i = 1;
for r=5:5:45
  Xpred_horpcac = tensor(zeros(sz));
  Xpred_l2e = tensor(zeros(sz));
  for j=1:10
    W = tensor(zeros(sz));
    W(cv.training(j))=1;
    result = wrapper_horpcac(X,W,r);
    Xpred_horpcac(cv.test(j))=result.X(cv.test(j));
    [T,tau] = tucker_l2e_opt(tensor(X),r,'init','tucker_als','W',W,'maxiters',5000,'taumax',50);
    T = tensor(T);
    Xpred_l2e(cv.test(j))=T(cv.test(j));
  end
  MADs_horpcac1(i) = sum(abs(double(Xpred_horpcac-X)),'all')/prod(sz);
  MADs_l2e1(i) = sum(abs(double(Xpred_l2e-X)),'all')/prod(sz);
  i = i+1;
end

reldiffs2_cpopt = zeros(9,1);
reldiffs2_tuckerals = zeros(9,1);
reldiffs2_horpcac = zeros(9,1);
reldiffs2_l2e = zeros(9,1);
W = tensor(ones(sz));
for j=1:1
  rng(j,'philox')
  %info = create_problem('Type','Tucker','Size',sz,'Num_Factors',[30 10 5],'Factor_Generator','orthogonal','Noise',0,'Core_Generator','randn');
  X = randn([50 50 50]);
  [truth,~,~,~,~,~] = hosvd(X,[30 10 5]);
  truth = tensor(truth);
  sz = size(X);
  Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
  stdtruth = std(truth(:));
  X = truth;
  X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
  for i=1:9
    r = i*5;
    Xhat_cpopt = cp_opt(X,r,'init','nvecs');
    reldiffs2_cpopt(i,j) = norm(tensor(Xhat_cpopt)-truth)/norm(truth);
    Xhat_tuckerals = tucker_als(X,r);
    reldiffs2_tuckerals(i,j) = norm(tensor(Xhat_tuckerals)-truth)/norm(truth);
    result_horpcac = wrapper_horpcac(X,'all_observed',r);
    Xhat_horpcac = result_horpcac.X;
    reldiffs2_horpcac(i,j) = norm(tensor(Xhat_horpcac)-truth)/norm(truth);
    Xhat_l2e = tucker_l2e_opt(X,r,'taumax',50);
    reldiffs2_l2e(i,j) = norm(tensor(Xhat_l2e)-truth)/norm(truth);
  end
end    

cv = cvpartition(prod(sz),'KFold',10);
MADs_horpcac2 = zeros(9,1);
MADs_l2e2 = zeros(9,1);
i = 1;
for r=5:5:45
  Xpred_horpcac = tensor(zeros(sz));
  Xpred_l2e = tensor(zeros(sz));
  for j=1:10
    W = tensor(zeros(sz));
    W(cv.training(j))=1;
    result = wrapper_horpcac(X,W,r);
    Xpred_horpcac(cv.test(j))=result.X(cv.test(j));
    [T,tau] = tucker_l2e_opt(tensor(X),r,'init','tucker_als','W',W,'maxiters',5000,'taumax',50);
    T = tensor(T);
    Xpred_l2e(cv.test(j))=T(cv.test(j));
  end
  MADs_horpcac2(i) = sum(abs(double(Xpred_horpcac-X)),'all')/prod(sz);
  MADs_l2e2(i) = sum(abs(double(Xpred_l2e-X)),'all')/prod(sz);
  i = i+1;
end

reldiffs3_cpopt = zeros(9,1);
reldiffs3_tuckerals = zeros(9,1);
reldiffs3_horpcac = zeros(9,1);
reldiffs3_l2e = zeros(9,1);
W = tensor(ones(sz));
for j=1:1
  rng(j,'philox')
  %info = create_problem('Type','Tucker','Size',sz,'Num_Factors',40,'Factor_Generator','orthogonal','Noise',0,'Core_Generator','randn');
  X = randn([50 50 50]);
  [truth,~,~,~,~,~] = hosvd(X,[35 35 35]);
  truth = tensor(truth);
  sz = size(X);
  Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
  stdtruth = std(truth(:));
  X = truth;
  X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
  for i=1:9
    r = i*5;
    Xhat_cpopt = cp_opt(X,r,'init','nvecs');
    reldiffs3_cpopt(i,j) = norm(tensor(Xhat_cpopt)-truth)/norm(truth);
    Xhat_tuckerals = tucker_als(X,r);
    reldiffs3_tuckerals(i,j) = norm(tensor(Xhat_tuckerals)-truth)/norm(truth);
    result_horpcac = wrapper_horpcac(X,'all_observed',r);
    Xhat_horpcac = result_horpcac.X;
    reldiffs3_horpcac(i,j) = norm(tensor(Xhat_horpcac)-truth)/norm(truth);
    Xhat_l2e = tucker_l2e_opt(X,r,'taumax',50);
    reldiffs3_l2e(i,j) = norm(tensor(Xhat_l2e)-truth)/norm(truth);
  end
end    

cv = cvpartition(prod(sz),'KFold',10);
MADs_horpcac3 = zeros(9,1);
MADs_l2e3 = zeros(9,1);
i = 1;
for r=5:5:45
  Xpred_horpcac = tensor(zeros(sz));
  Xpred_l2e = tensor(zeros(sz));
  for j=1:10
    W = tensor(zeros(sz));
    W(cv.training(j))=1;
    result = wrapper_horpcac(X,W,r);
    Xpred_horpcac(cv.test(j))=result.X(cv.test(j));
    [T,tau] = tucker_l2e_opt(tensor(X),r,'init','tucker_als','W',W,'maxiters',5000,'taumax',50);
    T = tensor(T);
    Xpred_l2e(cv.test(j))=T(cv.test(j));
  end
  MADs_horpcac3(i) = sum(abs(double(Xpred_horpcac-X)),'all')/prod(sz);
  MADs_l2e3(i) = sum(abs(double(Xpred_l2e-X)),'all')/prod(sz);
  i = i+1;
end

reldiffs1 = [mean(reldiffs1_cpopt,2) mean(reldiffs1_tuckerals,2) mean(reldiffs1_horpcac,2) mean(reldiffs1_l2e,2)]
reldiffs2 = [mean(reldiffs2_cpopt,2) mean(reldiffs2_tuckerals,2) mean(reldiffs2_horpcac,2) mean(reldiffs2_l2e,2)]
reldiffs3 = [mean(reldiffs3_cpopt,2) mean(reldiffs3_tuckerals,2) mean(reldiffs3_horpcac,2) mean(reldiffs3_l2e,2)]

csvwrite('reldiffs1_overfit.csv',reldiffs1)
csvwrite('reldiffs2_overfit.csv',reldiffs2)
csvwrite('reldiffs3_overfit.csv',reldiffs3)

a = [MADs_horpcac1 MADs_horpcac2 MADs_horpcac3]
b = [MADs_l2e1 MADs_l2e2 MADs_l2e3]
csvwrite('cv_horpca.csv',a)
csvwrite('cv_l2e.csv',b)


