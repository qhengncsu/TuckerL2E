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
addpath('./Tensor-Robust-Principal-Component-Analysis-TRPCA')
%note: all parameters are already tuned for optimal performance
%load in fMRI data (128*128*50)
load fMRI.mat;
rng(1,'philox');
sz = size(X);
Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
%add 25% of sparse noise
truth = double(X);
stdtruth = std(truth(:));
X(Oomega) = X(Oomega)+2*stdtruth*(rand(length(Oomega),1));

%the code for low-tubal-rank tRPCA (Liu et al. 2019)
opts.mu = 1e-4;
opts.tol = 1e-5;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;
[n1,n2,n3] = size(double(X));
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat_trpca,E,err,iter] = trpca_tnn(double(X),lambda,opts);

%print out the relative error for tRPCA
reldiff_trpca = norm(tensor(Xhat_trpca)-tensor(truth))/norm(tensor(truth))

%the code for HoRPCA-C (Goldfarb and Qin, 2014)
result_horpcac = wrapper_horpcac(tensor(X),'all_observed',[96 96 15]);
Xhat_horpcac = result_horpcac.X;

%print out the relative error for HoRPCA-C
reldiff_horpcac = norm(tensor(Xhat_horpcac)-tensor(truth))/norm(tensor(truth))

%the code for RGrad (Cai et al. 2022)
opts.mu0 = 1;
[Xhat_rgrad,S,time] = tRPCA_RGrad(double(X),[96 96 15],1,0.27,double(X),opts);

%print out the relative error for RGrad
reldiff_rgrad = norm(tensor(Xhat_rgrad)-tensor(truth))/norm(tensor(truth))

%the code for Tucker-L2E (our method) (warning: this will take about 3-5 minutes)
T = tucker_l2e_opt(tensor(X),[128 128 20],'maxiters',5000,'taumax',50);
Xhat_l2e = tensor(T);

%print out the relative error for Tucker-L2E
reldiff_l2e = norm(tensor(T)-tensor(truth))/norm(tensor(truth))

%what happens if we use rank (128,128,20) for HoRPCA-C and RGrad?
%result_horpcac = wrapper_horpcac(tensor(X),'all_observed',[128 128 20]);
%Xhat_horpcac_128 = result_horpcac.X;
%reldiff_horpcac_128 = norm(tensor(Xhat_horpcac_128)-tensor(truth))/norm(tensor(truth))

%[Xhat_rgrad_128,S,time] = tRPCA_RGrad(double(X),[128 128 20],1,0.27,double(X),opts);
%reldiff_rgrad_128 = norm(tensor(Xhat_rgrad_128)-tensor(truth))/norm(tensor(truth))
%they are much less accurate due to overfitting to the sparse noise, in 
%in other words, they are not stable at rank (128,128,20)

%visualize the recovery results
reldiff_noisy = norm(tensor(X)-tensor(truth))/norm(tensor(truth))
set(gcf,'Position',[100 100 1200 800])
subplot(2,3,1)
imshow(uint16(double(truth(:,:,30))),'DisplayRange',[]);
title('Original','FontSize',15)
subplot(2,3,2)
imshow(uint16(double(X(:,:,30))),'DisplayRange',[]);
title(strcat('Noisy, RE=',sprintf('%.4f',reldiff_noisy)),'FontSize',15)
subplot(2,3,3)
imshow(uint16(double(Xhat_trpca(:,:,30))),'DisplayRange',[]);
title(strcat('tRPCA, RE=',sprintf('%.4f',reldiff_trpca)),'FontSize',15)
subplot(2,3,4)
imshow(uint16(double(Xhat_horpcac(:,:,30))),'DisplayRange',[]);
title(strcat('HoRPCA-C, RE=',sprintf('%.4f',reldiff_horpcac)),'FontSize',15)
subplot(2,3,5)
imshow(uint16(double(Xhat_rgrad(:,:,30))),'DisplayRange',[]);
title(strcat('RGrad, RE=',sprintf('%.4f',reldiff_rgrad)),'FontSize',15)
subplot(2,3,6)
imshow(uint16(double(Xhat_l2e(:,:,30))),'DisplayRange',[]);
title(strcat('Tucker-L2E, RE=',sprintf('%.4f',reldiff_l2e)),'FontSize',15)