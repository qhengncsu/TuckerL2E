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

%load data, please ignore the warnings
load dorrit
%remove certain regions of the dataset
full = 1:27;
outlying_samples=[2 3 5 10];
partial = setdiff(full,outlying_samples);
X = X1(partial,:,7:24);

%code for RGrad
opts.mu0 = 5;
[Xhat_rgrad,S,time] = tRPCA_RGrad(double(X),[4 4 4],1,0.11,double(X),opts);

%code for HoRPCA-C (Goldfarb and Qin, 2014)
result_horpcac = wrapper_horpcac(tensor(X),'all_observed',4);
Xhat_horpcac = result_horpcac.X;

%code for Tucker-L2E (our method)
Xhat_l2e = tucker_l2e_opt(tensor(X),4,'init','tucker_als','taumax',50,'maxiters',5000);
Xhat_l2e = tensor(Xhat_l2e);

%use CP-OPT to extract the CP factors
M1 = cp_opt(tensor(X),4,'init','nvecs');
M2 = cp_opt(tensor(Xhat_rgrad),4,'init','nvecs');
M3 = cp_opt(Xhat_horpcac,4,'init','nvecs');
M4 = cp_opt(Xhat_l2e,4,'init','nvecs');

%visualize the recovery results
emission_scale = 241:2:481;
excitation_scale = 230:5:315;
emission1 = M1.U{2};
excitation1 = M1.U{3};
emission2 = M2.U{2};
excitation2 = M2.U{3};
emission3 = M3.U{2};
excitation3 = M3.U{3};
emission4 = M4.U{2};
excitation4 = M4.U{3};
t=tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[100 100 1200 600])
nexttile
plot(emission_scale,emission1(:,1),'LineWidth',2,'DisplayName','tryptophan')
hold on 
plot(emission_scale,emission1(:,3),'LineWidth',2,'DisplayName','dopa')
hold on 
plot(emission_scale,emission1(:,2),'LineWidth',2,'DisplayName','hydroquinone')
hold on 
plot(emission_scale,emission1(:,4),'LineWidth',2,'DisplayName','phenylalanine')
hold off
xlim([241 481])
title('Noisy emission loadings')
l = legend('show','Location','northeast')
nexttile
plot(emission_scale,emission2(:,1),'LineWidth',2)
hold on 
plot(emission_scale,emission2(:,3),'LineWidth',2)
hold on 
plot(emission_scale,emission2(:,2),'LineWidth',2)
hold on 
plot(emission_scale,emission2(:,4),'LineWidth',2)
hold off
xlim([241 481])
title('RGrad emission loadings')

nexttile
plot(emission_scale,emission3(:,1),'LineWidth',2)
hold on 
plot(emission_scale,emission3(:,2),'LineWidth',2)
hold on 
plot(emission_scale,emission3(:,3),'LineWidth',2)
hold on 
plot(emission_scale,emission3(:,4),'LineWidth',2)
hold off
xlim([241 481])
title('HoRPCA-C emission loadings')

nexttile
plot(emission_scale,emission4(:,1),'LineWidth',2)
hold on 
plot(emission_scale,emission4(:,2),'LineWidth',2)
hold on 
plot(emission_scale,emission4(:,3),'LineWidth',2)
hold on 
plot(emission_scale,emission4(:,4),'LineWidth',2)
hold off
xlim([241 481])
title('Tucker-L2E emission loadings')

nexttile
plot(excitation_scale,excitation1(:,1),'LineWidth',2')
hold on 
plot(excitation_scale,excitation1(:,3),'LineWidth',2)
hold on 
plot(excitation_scale,excitation1(:,2),'LineWidth',2)
hold on 
plot(excitation_scale,excitation1(:,4),'LineWidth',2)
hold off
xlim([230 315])
title('Noisy excitation loadings')

nexttile
plot(excitation_scale,excitation2(:,1),'LineWidth',2)
hold on 
plot(excitation_scale,excitation2(:,3),'LineWidth',2)
hold on 
plot(excitation_scale,excitation2(:,2),'LineWidth',2)
hold on 
plot(excitation_scale,excitation2(:,4),'LineWidth',2)
hold off
xlim([230 315])
ylim([-0.02 0.6])
title('RGrad excitation loadings')

nexttile
plot(excitation_scale,excitation3(:,1),'LineWidth',2)
hold on 
plot(excitation_scale,excitation3(:,2),'LineWidth',2)
hold on 
plot(excitation_scale,excitation3(:,3),'LineWidth',2)
hold on 
plot(excitation_scale,excitation3(:,4),'LineWidth',2)
hold off
xlim([230 315])
%ylim([-0.02 0.6])
title('HoRPCA-C excitation loadings')

nexttile
plot(excitation_scale,excitation4(:,1),'LineWidth',2)
hold on 
plot(excitation_scale,excitation4(:,2),'LineWidth',2)
hold on 
plot(excitation_scale,excitation4(:,3),'LineWidth',2)
hold on 
plot(excitation_scale,excitation4(:,4),'LineWidth',2)
hold off
xlim([230 315])
title('Tucker-L2E excitation loadings')
xlabel(t,'Wavelength');
ylabel(t,'Intensity');

