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
reldiffs_tals = zeros(9,50);
reldiffs_rpca = zeros(9,50);
reldiffs_horpcas = zeros(9,50);
reldiffs_horpcac = zeros(9,50);
reldiffs_rgrad = zeros(9,50);
reldiffs_l2e = zeros(9,50);
for i=1:9
  for j=1:50
    rng(i*50+j,'philox')
    %info = create_problem('Type', 'Tucker', 'Size',[50 50 50],'Num_Factors',5*i,'Noise', 0.1,'M',0.0,'Factor_Generator','rand','Core_Generator','randn');
    %X = info.Data;
    X = randn([50 50 50]);
    [truth,~,~,~,~,~] = hosvd(X,[5*i 5*i 5*i]);
    truth = tensor(truth);
    sz = size(X);
    Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
    stdtruth = std(truth(:));
    X = truth;
    X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
    E = tensor(randn(sz)*0.0*stdtruth);
    X = X+E;
    opts.mu0 = 5;
    [Xhat_rgrad,~,time] = tRPCA_RGrad(double(X),[5*i 5*i 5*i],1,0.25,double(X),opts);
    reldiffs_rgrad5 = norm(tensor(Xhat_rgrad)-truth)/norm(truth);
    opts.mu0 = 1;
    [Xhat_rgrad,~,time] = tRPCA_RGrad(double(X),[5*i 5*i 5*i],1,0.25,double(X),opts);
    reldiffs_rgrad1 = norm(tensor(Xhat_rgrad)-truth)/norm(truth);
    reldiffs_rgrad(i,j) = min([reldiffs_rgrad5 reldiffs_rgrad1]);
    [Xhat_horpcas min_reldiff] = twg_horpcas(X,'all_observed',truth); 
    reldiffs_horpcas(i,j) = min_reldiff;
    result_horpcac = wrapper_horpcac(X,'all_observed',5*i);
    Xhat_horpcac = result_horpcac.X;
    reldiffs_horpcac(i,j) = norm(truth-tensor(Xhat_horpcac))/norm(truth);
    [T,tau] = tucker_l2e_opt(tensor(X),5*i,'taumax',50);
    reldiffs_l2e(i,j) = norm(tensor(T)-truth)/norm(truth);
    %Xhat_cpopt = cp_opt(X,5*i,'init','nvecs');
    %reldiffs_cpopt(i,j) = norm(truth-Xhat_cpopt)/norm(truth);
    Xhat_tals = tucker_als(X,5*i);
    reldiffs_tals(i,j) = norm(truth-tensor(Xhat_tals))/norm(truth);
    %[model] = BayesRCP(X, 'init', 'ml', 'initVar', 1, 'maxRank', 100, 'dimRed', 1);
    %Xhat_brtf = tensor(ktensor(model.Z));
    %reldiffs_brtf(i,j) = norm(truth-tensor(Xhat_brtf))/norm(truth);
    [Xhat_rpca min_reldiff] = twg_rpca(X,'all_observed',truth); 
    reldiffs_rpca(i,j) = min_reldiff;
  end
end
result = [mean(reldiffs_tals,2) mean(reldiffs_rpca,2) mean(reldiffs_horpcas,2) mean(reldiffs_horpcac,2) mean(reldiffs_rgrad,2) mean(reldiffs_l2e,2)]
csvwrite('result_tucker_gross_nodense.csv',result)