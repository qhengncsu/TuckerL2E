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
reldiffs_horpca = zeros(11,11,20);
reldiffs_l2e = zeros(11,11,20);
for k=1:20
  for i=1:11
    rng(i*10+k,'multFibonacci')
    %info = create_problem('Type', 'Tucker', 'Size',[50 50 50],'Num_Factors',25+2*(i-1),'Noise', 0.0,'M',0.0,'Factor_Generator','orthogonal','Core_Generator','randn');
    X = randn([50 50 50]);
    [truth,~,~,~,~,~] = hosvd(X,[25+2*(i-1) 25+2*(i-1) 25+2*(i-1)]);
    truth = tensor(truth);
    sz = size(X);
    stdtruth = std(truth(:));
    for j=1:11
      rng(i*11+j,'multFibonacci')
      Oomega = randsample(prod(sz), int64(round(prod(sz)*0.05*(j-1))));
      X = truth;
      X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
      result = wrapper_horpcac(X,'all_observed',25+2*(i-1));
      reldiffs_horpca(i,j,k) = norm(tensor(result.X)-truth)/norm(truth);
      [T,tau] = tucker_l2e_opt(tensor(X),25+2*(i-1),'taumax',100,'maxiters',5000);
      reldiffs_l2e1 = norm(tensor(T)-truth)/norm(truth);
      [T,tau] = tucker_l2e_opt(tensor(X),25+2*(i-1),'taumax',50,'maxiters',5000);
      reldiffs_l2e2 = norm(tensor(T)-truth)/norm(truth);
      reldiffs_l2e(i,j,k) = min([reldiffs_l2e1 reldiffs_l2e2]);
    end
  end
  a = mean(reldiffs_horpca(:,:,1:k),3)
  b = mean(reldiffs_l2e(:,:,1:k),3)
  csvwrite('pt_horpca_10std_tucker.csv',a)
  csvwrite('pt_l2e_10std_tucker.csv',b)
end