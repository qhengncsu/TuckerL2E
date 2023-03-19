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
    info = create_problem('Type', 'CP', 'Size',[50 50 50],'Num_Factors',25+2*(i-1),'Noise', 0.0,'M',0,'Factor_Generator','randn');
    for j=1:11
      X = info.Data;
      sz = size(X);
      truth = tensor(info.Soln);
      Oomega = randsample(prod(sz), int64(round(prod(sz)*0.05*(j-1))));
      stdtruth = std(truth(:));
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
  a = mean(reldiffs_horpca(:,:,1:k),3);
  b = mean(reldiffs_l2e(:,:,1:k),3);
  csvwrite('pt_horpca_10std_cp.csv',a)
  csvwrite('pt_l2e_10std_cp.csv',b)
end
