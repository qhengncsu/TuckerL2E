load('it4tnsr_notime');     %%%%% load the data
A = log(1+it4tnsr_notime);   %%%% log transformation
load('country_region')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameter setting
maxit = 70;
alpha = 0.3;  %% S figures are obtained with alpha=0.03
gamma = 1;
beta = 1;
mu0 = 5;
r =[3,3,3];
stop_thres = 0.0001;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   initialization by HOSVD
init_hosvd = 1;
if init_hosvd == 0
    perturb_sigma = argin.perturb_sigma;
    %rng(2);
    perturb = randn(d)*perturb_sigma;
    T = Ttrue + perturb; 
    [T,C,U,UT,~,~] = hosvd(T,r); 
else
    [T,C,U,UT,~,~] = hosvd(A,r); 
end  
S = threshold(T,gamma,alpha,A); %%%% initialization for S

init_err = norm(A(:) - T(:)-S(:));
disp(['initialization rel err.: ', num2str(init_err)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% iteration of RGrad algorithm
tic;
Tprev = T;
for i = 1:maxit
    G = T + S - A;
    [G_p,~,~] = mani_proj(G,C,U,UT);
    W = T - beta*G_p;
    [T,C,U,UT,~,~] = hosvd(W,r); 
    tildeT = trimtensor(T,mu0); % norm(tildeT(:))
    S = threshold(tildeT,gamma,alpha,A); 
    
    relerrList(i,1) = norm(A(:) - T(:)-S(:)); 
    
    disp(['iter:', num2str(i),', obj val:',num2str(relerrList(i,1))]);
    if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres
        break;
    end
    Tprev = tildeT;
end
runtime = toc




%%%%% MDS and plot of exports.  Image size, width=833 and height=745, pos1=1000, pos2=593
% pU=cmdscale(squareform(pdist(U{1})),2);
% figure;
% africa_id=find(country_region=='Africa');
% asia_id=find(country_region=='Asia');
% america_id=find(country_region=='America');
% europe_id=find(country_region=='Europe');
% plot(pU(:,1),pU(:,2),'.','MarkerSize',14)
% text(pU(europe_id,1),pU(europe_id,2),country_names(europe_id),'color','b', 'fontsize',14)
% text(pU(america_id,1),pU(america_id,2),country_names(america_id),'color','r','fontsize',14)
% text(pU(asia_id,1),pU(asia_id,2),country_names(asia_id),'color','m','fontsize',14)
% text(pU(africa_id,1),pU(africa_id,2),country_names(africa_id),'color','k','fontsize',14)
% xlabel("First Component after Multidimensional Scaling",'fontsize',14);
% ylabel("Second Component after Multidimensional Scaling",'fontsize',14);



%%%%% MDS and plot of imports.  Image size, width=833 and height=745
pU=cmdscale(squareform(pdist(U{2})),2);
figure;
africa_id=find(country_region=='Africa');
asia_id=find(country_region=='Asia');
america_id=find(country_region=='America');
europe_id=find(country_region=='Europe');
%tmp=pU(:,2);pU(:,2)=pU(:,1);pU(:,1)=tmp;
plot(pU(:,1),pU(:,2),'.','MarkerSize',14)
set(gcf,'Position',[1000 593 833 745])
text(pU(europe_id,1),pU(europe_id,2),country_names(europe_id),'color','b', 'fontsize',14)
text(pU(america_id,1),pU(america_id,2),country_names(america_id),'color','r','fontsize',14)
text(pU(asia_id,1),pU(asia_id,2),country_names(asia_id),'color','m','fontsize',14)
text(pU(africa_id,1),pU(africa_id,2),country_names(africa_id),'color','k','fontsize',14)
xlabel("First Component after Multidimensional Scaling",'fontsize',14);
ylabel("Second Component after Multidimensional Scaling",'fontsize',14);


%%%%% MDS and plot of commodities.  Image size, width=833 and height=745
% pU=cmdscale(squareform(pdist(U{3})),2);
% figure;
% plot(pU(:,1),pU(:,2),'.','MarkerSize',14)
% set(gcf,'Position',[1000 593 833 745])
% text(pU(:,1),pU(:,2),num2str(1:100),'color','b', 'fontsize',14)
% xlabel("First Component after Multidimensional Scaling",'fontsize',14);
% ylabel("Second Component after Multidimensional Scaling",'fontsize',14);


% %%%%% Heatmaps of S
% pS=abs(S(:,:,2));
% pS=pS-diag(diag(pS));
% figure
% hS1=heatmap(pS);
% hS1.XDisplayLabels = country_names;
% hS1.YDisplayLabels = country_names;
% 
% 
% 
% %%%%% Heatmaps of S
% pS=abs(S(:,:,27));
% pS=pS-diag(diag(pS));
% figure
% hS2=heatmap(pS);
% hS2.XDisplayLabels = country_names;
% hS2.YDisplayLabels = country_names;
% 
% 
% 
% %%%%% Heatmaps of S
% pS=abs(S(:,:,31));
% pS=pS-diag(diag(pS));
% figure
% hS3=heatmap(pS);
% hS3.XDisplayLabels = country_names;
% hS3.YDisplayLabels = country_names;
% 
% 
% 
% %%%%% Heatmaps of S
% pS=abs(S(:,:,24));
% pS=pS-diag(diag(pS));
% figure
% hS4=heatmap(pS);
% hS4.XDisplayLabels = country_names;
% hS4.YDisplayLabels = country_names;
% 
% 
% for i=71:100
% pS=abs(S(:,:,i));
% pS=pS-diag(diag(pS));
% figure
% hS1=heatmap(pS);
% hS1.XDisplayLabels = country_names;
% hS1.YDisplayLabels = country_names;
% title(['slice: ',num2str(i)])
% end

%% import S slices 23, 27, 17, 13, 55, 52, 45, 80, 77, 99, 98, 65, 47, 96
%pos=[271   398   953   811]
%sind=[23, 27, 17, 13, 55, 52, 45, 80, 77, 99, 98, 65, 47, 96];
%product_name=["Food Industries","Mineral Fuels","Sugar","LAC, GUMS, RESINS","MAN-MADE STAPLE FIBRES","Cotton","Cork","Tin","77 not found","SPECIAL PURPOSE","Miscellaneous Not Classified","Head Gear","Pulp of Wood","MISCELLANEOUS MANUFACTURED ARTICLES"];
%for i=1:length(sind)
%pS=abs(S(:,:,sind(i)));
%pS=pS-diag(diag(pS));
%figure
%hS4=heatmap(pS,'CellLabelColor','none');
%hS4.XDisplayLabels = country_names;
%hS4.YDisplayLabels = country_names;
%title(strcat("Product ",num2str(sind(i)),": ",product_name(i)))
%set(gcf,'Position',[271   398   953   811])
%end


%% Slices to draw 98, 27, 55, 99, 80, 45, 52, 17

