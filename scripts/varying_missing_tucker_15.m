addpath('./HoRPCA/inexact_alm_rpca')
addpath('./HoRPCA/data')
addpath('./HoRPCA/lightspeed')
addpath('./HoRPCA/PROPACK')
addpath('./HoRPCA/code/rpca')
addpath('./HoRPCA/code/tc')
addpath('./HoRPCA/code/utils')
addpath('./HoRPCA')
addpath('./HoRPCA/code/utils')
addpath('./HoRPCA')
addpath('./BRTF/BRTF')
addpath('./BRTF/BRTF/mykhatrirao')
addpath('./BRTF/BRTF/tensor_plot')
addpath('./TuckerL2E')
addpath('./TuckerL2E/L-BFGS-B-C/Matlab')
addpath('./TuckerL2E/tensor_toolbox-v3.2.1')
addpath('./RGrad')
addpath('./RGrad/functions')
reldiffs_cpopt = zeros(9,20);
reldiffs_tals = zeros(9,20);
reldiffs_brtf = zeros(9,20);
reldiffs_rpca = zeros(9,20);
reldiffs_horpcas = zeros(9,20);
reldiffs_horpcac = zeros(9,20);
reldiffs_rgrad = zeros(9,20);
reldiffs_l2e_spectral = zeros(9,20);
reldiffs_l2e_tuckerals = zeros(9,20);
for j=1:20
  for i=1:18
    rng(i*20+j,'multFibonacci')
    X = randn([50 50 50]);
    [truth,~,~,~,~,~] = hosvd(X,[15 15 15]);
    truth = tensor(truth);
    sz = size(X);
    Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
    stdtruth = std(truth(:));
    X = truth;
    X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
    E = tensor(randn(sz)*0.1*stdtruth);
    X = X+E;
    W = tensor(ones(sz));
    Womega = randsample(prod(sz), int64(round(prod(sz)*(i*0.05))));
    W(Womega) = 0;
    [Xhat_horpcas min_reldiff] = twg_horpcas(X,W,truth); 
    reldiffs_horpcas(i,j) = min_reldiff;
    result_horpcac = wrapper_horpcac(X,W,15);
    Xhat_horpcac = result_horpcac.X;
    reldiffs_horpcac(i,j) = norm(truth-tensor(Xhat_horpcac))/norm(truth);
    [T,tau] = tucker_l2e_opt(X,15,'taumax',50,'W',W,'init','spectral');
    reldiffs_l2e_spectral(i,j) = norm(tensor(T)-truth)/norm(truth);
    [T,tau] = tucker_l2e_opt(X,15,'taumax',50,'W',W);
    reldiffs_l2e_tuckerals(i,j) = norm(tensor(T)-truth)/norm(truth);
    %[model] = BayesRCP(double(X), 'init', 'ml', 'initVar', 1, 'maxRank', 100, 'dimRed', 1);
    %Xhat_brtf = tensor(ktensor(model.Z));
    %reldiffs_brtf(i,j) = norm(truth-tensor(Xhat_brtf))/norm(truth);
    %[Xhat_rpca min_reldiff] = twg_rpca(X,'all_observed',truth); 
    %reldiffs_rpca(i,j) = min_reldiff;
  end
  result = [mean(reldiffs_horpcas(:,1:j),2) mean(reldiffs_horpcac(:,1:j),2) mean(reldiffs_l2e_spectral(:,1:j),2) mean(reldiffs_l2e_tuckerals(:,1:j),2)]
  csvwrite('result_missing_tucker_15.csv',result)
end