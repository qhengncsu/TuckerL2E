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

rng(1,'philox')
reldiffs_horpcac = zeros(9,1);
reldiffs_rgrad = zeros(9,1);
reldiffs_l2e = zeros(9,1);
for i=1:9
    %generate a tensor with Tucker rank (5i,5i,5i) (i=1,2,3,...,9)
    X = randn([50 50 50]);
    [truth,~,~,~,~,~] = hosvd(X,[5*i 5*i 5*i]);
    truth = tensor(truth);

    %add 25% of large outliers
    sz = size(X);
    Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
    stdtruth = std(truth(:));
    X = truth;
    X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);

    %add dense normal noise of relative scale 0.1
    E = tensor(randn(sz)*0.1*stdtruth);
    X = X+E;

    %code for HoRPCA-C (Goldfarb and Qin, 2014)
    result_horpcac = wrapper_horpcac(X,'all_observed',[5*i 5*i 5*i]);
    Xhat_horpcac = result_horpcac.X;

    %record the relative error for HoRPCA-C
    reldiffs_horpcac(i) = norm(truth-tensor(Xhat_horpcac))/norm(truth);

    %code for RGrad (Cai et al. 2022)
    opts.mu0 = 5;
    [Xhat_rgrad,~,time] = tRPCA_RGrad(double(X),[5*i 5*i 5*i],1,0.25,double(X),opts);

    %print out the relative error for RGrad
    reldiffs_rgrad(i) = norm(tensor(Xhat_rgrad)-truth)/norm(truth);

    %code for Tucker-L2E (our method)
    [T,tau] = tucker_l2e_opt(tensor(X),[5*i 5*i 5*i],'taumax',50);

    %print out the relative error for Tucker-L2E
    reldiffs_l2e(i) = norm(tensor(T)-truth)/norm(truth);
end

%visualize the recovery results
ranks = 5:5:45;
t=tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[100 100 640 640])
nexttile
line1 = plot(ranks,reldiffs_horpcac,'-sr','DisplayName','HoRPCA-C','LineWidth',1.5);
hold on 
line2 = plot(ranks,reldiffs_rgrad,'-x','DisplayName','RGrad','color',[0.5 0 0.8],'LineWidth',1.5);
hold on
line3 = plot(ranks,reldiffs_l2e,'-ob','DisplayName','Tucker-L2E','LineWidth',1.5);
hold off
l = legend('show','Location','northwest')
xlim([5 45])
xlabel(t,'Rank');
ylabel(t,'Relative Error');