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
addpath('./BRTF/BRTF')
addpath('./BRTF/BRTF/mykhatrirao')
addpath('./BRTF/BRTF/tensor_plot')

N = 1;
reldiffs_cpopt = zeros(9,N);
reldiffs_brtf = zeros(9,N);
reldiffs_horpcac = zeros(9,N);
reldiffs_rgrad = zeros(9,N);
reldiffs_l2e = zeros(9,N);
for j=1:N
    rng('philox',j)
    for i=1:9
        %generate a tensor with CP rank (5i,5i,5i) (i=1,2,3,...,9)
        info = create_problem('Type', 'CP', 'Size',[50 50 50],'Num_Factors',5*i,'Noise', 0.0,'M',0.0,'Factor_Generator','randn');
        X = info.Data;
        sz = size(X);
        truth = tensor(info.Soln);
    
        %add 25% of large outliers
        Oomega = randsample(prod(sz), int64(round(prod(sz)*0.25)));
        stdtruth = std(truth(:));
        X = truth;
        X(Oomega) = X(Oomega)+10*stdtruth*(2*rand(length(Oomega),1)-1);
    
        %add dense normal noise of relative scale 0.1
        E = tensor(randn(sz)*0.1*stdtruth);
        X = X+E;
        
        %code for CP-OPT
        Xhat_cpopt =cp_opt(X,5*i,'init','nvecs');
        
        %record the relative error for Tucker-ALS
        reldiffs_cpopt(i,j) = norm(truth-tensor(Xhat_cpopt))/norm(truth);
        
        %code for BRTF (Zhao et al., 2014)
        [model] = BayesRCP(double(X), 'init', 'ml', 'initVar', 1, 'maxRank', 100, 'dimRed', 1);
        Xhat_brtf = tensor(ktensor(model.Z));
        
        %record the relative error for HoRPCA-S
        reldiffs_brtf(i,j) = norm(truth-tensor(Xhat_brtf))/norm(truth);
    
        %code for HoRPCA-C (Goldfarb and Qin, 2014)
        result_horpcac = wrapper_horpcac(X,'all_observed',[5*i 5*i 5*i]);
        Xhat_horpcac = result_horpcac.X;
    
        %record the relative error for HoRPCA-C
        reldiffs_horpcac(i,j) = norm(truth-tensor(Xhat_horpcac))/norm(truth);
    
        %code for RGrad (Cai et al. 2022)
        opts.mu0 = 5;
        [Xhat_rgrad,~,time] = tRPCA_RGrad(double(X),[5*i 5*i 5*i],1,0.25,double(X),opts);
    
        %print out the relative error for RGrad
        reldiffs_rgrad(i,j) = norm(tensor(Xhat_rgrad)-truth)/norm(truth);
    
        %code for Tucker-L2E (our method)
        [T,tau] = tucker_l2e_opt(tensor(X),[5*i 5*i 5*i],'taumax',50);
    
        %print out the relative error for Tucker-L2E
        reldiffs_l2e(i,j) = norm(tensor(T)-truth)/norm(truth);
    end
end
%cap the relative errors at 1
reldiffs_cpopt(reldiffs_cpopt>1) = 1;
reldiffs_brtf(reldiffs_brtf>1) = 1;
reldiffs_horpcac(reldiffs_horpcac>1) = 1;
reldiffs_rgrad(reldiffs_rgrad>1) = 1;
reldiffs_l2e(reldiffs_l2e>1) = 1;

%visualize the recovery results
ranks = 5:5:45;
set(gcf,'Position',[100 100 800 600])
line1 = plot(ranks,mean(reldiffs_cpopt,2),'-+','DisplayName','CP-OPT','color','#77AC30','LineWidth',1.5);
hold on
line2 = plot(ranks,mean(reldiffs_brtf,2),'-d','DisplayName','BRTF','color','#4DBEEE','LineWidth',1.5);
hold on
line3 = plot(ranks,mean(reldiffs_horpcac,2),'-sr','DisplayName','HoRPCA-C','LineWidth',1.5);
hold on 
line4 = plot(ranks,mean(reldiffs_rgrad,2),'-x','DisplayName','RGrad','color',[0.5 0 0.8],'LineWidth',1.5);
hold on
line5 = plot(ranks,mean(reldiffs_l2e,2),'-ob','DisplayName','Tucker-L2E','LineWidth',1.5);
hold off
l = legend('show','Location','northwest');
xlim([5 45])
xlabel('Rank');
ylabel('Relative Error');