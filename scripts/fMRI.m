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
addpath('./Tensor-Robust-Principal-Component-Analysis-TRPCA-master')

load fMRI.mat;
rng(1,'philox');
sz = size(X);
Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
truth = double(X);
stdtruth = std(truth(:));
X(Oomega) = X(Oomega)+2*stdtruth*(rand(length(Oomega),1));
sz = size(X);
holdout = randsample(prod(sz), int64(round(prod(sz)*0.1)));
W = tensor(ones(sz));
W(holdout)=0;
T = tucker_l2e_opt(tensor(X),[64 64 10],'W',W,'maxiters',5000,'taumax',50);
T = tensor(T);
sum(abs(T(holdout)-X(holdout)),'all')/length(holdout)

T = tucker_l2e_opt(tensor(X),[96 96 15],'W',W,'maxiters',5000,'taumax',50);
T = tensor(T);
sum(abs(T(holdout)-X(holdout)),'all')/length(holdout)

T = tucker_l2e_opt(tensor(X),[128 128 20],'W',W,'maxiters',5000,'taumax',50);
T = tensor(T);
sum(abs(T(holdout)-X(holdout)),'all')/length(holdout)

result_horpcac = wrapper_horpcac(tensor(X),W,[96 96 15]);
Xhat_horpcac = result_horpcac.X;
sum(abs(Xhat_horpcac(holdout)-X(holdout)),'all')/length(holdout)

result_horpcac = wrapper_horpcac(tensor(X),W,[128 128 20]);
Xhat_horpcac = result_horpcac.X;
sum(abs(Xhat_horpcac(holdout)-X(holdout)),'all')/length(holdout)

result_horpcac = wrapper_horpcac(tensor(X),W,[64 64 10]);
Xhat_horpcac = result_horpcac.X;
sum(abs(Xhat_horpcac(holdout)-X(holdout)),'all')/length(holdout)

bics = zeros(11,1);
alphas = 0.2:0.01:0.3;
opts.mu0 = 1;
for i=1:11
    [Xhat_rgrad,S,time] = tRPCA_RGrad(double(X),[96 96 15],1,alphas(i),double(X),opts);
    bics(i) = bic(double(X),S,Xhat_rgrad,[96 96 15]);
end
alpha = alphas(argmin(bics));
[Xhat_rgrad,S,time] = tRPCA_RGrad(double(X),[96 96 15],1,alpha,double(X),opts);
reldiff_rgrad = norm(tensor(Xhat_rgrad)-tensor(truth))/norm(tensor(truth))

result_horpcac = wrapper_horpcac(tensor(X),'all_observed',[96 96 15]);
Xhat_horpcac = result_horpcac.X;
reldiff_horpcac = norm(tensor(Xhat_horpcac)-tensor(truth))/norm(tensor(truth))

T = tucker_l2e_opt(tensor(X),[128 128 20],'maxiters',5000,'taumax',50);
reldiff_l2e = norm(tensor(T)-tensor(truth))/norm(tensor(truth))
Xhat_l2e = tensor(T);

opts.mu = 1e-4;
opts.tol = 1e-5;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;
[n1,n2,n3] = size(double(X));
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat_trpca,E,err,iter] = trpca_tnn(double(X),lambda,opts);
reldiff_trpca = norm(tensor(Xhat_trpca)-tensor(truth))/norm(tensor(truth))

reldiff_noisy = norm(tensor(X)-tensor(truth))/norm(tensor(truth))
t=tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[0 0 1200 960])
nexttile
imshow(uint16(double(truth(:,:,30))),'DisplayRange',[]);
title('Original','FontSize',20)
nexttile
imshow(uint16(double(X(:,:,30))),'DisplayRange',[]);
title(strcat('Noisy, RE=',sprintf('%.4f',reldiff_noisy)),'FontSize',20)
nexttile
imshow(uint16(double(Xhat_trpca(:,:,30))),'DisplayRange',[]);
title(strcat('tRPCA, RE=',sprintf('%.4f',reldiff_trpca)),'FontSize',20)
nexttile
imshow(uint16(double(Xhat_horpcac(:,:,30))),'DisplayRange',[]);
title(strcat('HoRPCA-C, RE=',sprintf('%.4f',reldiff_horpcac)),'FontSize',20)
nexttile
imshow(uint16(double(Xhat_rgrad(:,:,30))),'DisplayRange',[]);
title(strcat('RGrad, RE=',sprintf('%.4f',reldiff_rgrad)),'FontSize',20)
nexttile
imshow(uint16(double(Xhat_l2e(:,:,30))),'DisplayRange',[]);
title(strcat('Tucker-L2E, RE=',sprintf('%.4f',reldiff_l2e)),'FontSize',20)
