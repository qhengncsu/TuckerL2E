load('Adj_nodummy.mat')     %%%%%% load the dataset
load('author')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter setting
alpha = 0.0005;  
gamma = 1;
beta = 0.8;
mu0 = 5;
maxit = 10;
k=4;
r =[k,k,k]; 
stop_thres = 0.001;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization by HOSVD
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
S = threshold(T,gamma,alpha,A); 
init_err = norm(A(:) - T(:)-S(:));
disp(['initialization rel err.: ', num2str(init_err)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% iteration of RGrad algorithm
tic;
Tprev = T;
for i = 1:maxit
    G = T + S - A;
    [G_p,~,~] = mani_proj(G,C,U,UT);
    W = T - beta*G_p;
    [T,C,U,UT,~,~] = hosvd(W,r); 
    tildeT = trimtensor(T,mu0);
    S = threshold(tildeT,gamma,alpha,A); 
    relerrList(i,1) = norm(A(:) - T(:)-S(:)); 
    disp(['iter:', num2str(i),', obj val:',num2str(relerrList(i,1))]);
    if norm(T(:) - Tprev(:))/norm(T(:)) < stop_thres
        break;
    end
    Tprev = tildeT;
end
runtime = toc




%%%%% SCORE and MDS and plot of authors relation -- alpha=0.0001
tmpU=U{1};
tmp=tmpU(:,2:k)./tmpU(:,1);tmp(isnan(tmp))=100;delta=15;tmp=trim(tmp,delta);
[cidx, ctrs] = kmeans(tmp, 3, 'Replicates',20);
g1=find(cidx==1);
g2=find(cidx==2);
g3=find(cidx==3);

%%%%%%%% too many authors, just display some of them
key_author=["David Dunson","Amy H Herring","Anirban Bhattacharya","Lawrence Carin","Yongtao Guan","Brian J Reich","Song Xi Chen","Fang Yao","Ying Wei","Heping Zhang","Harrison H Zhou","Debashis Paul","Yi Lin","Michael J Todd","Kani Chen","Qiwei Yao","Qi-Man Shao","Ming Yuan","Richard Samworth","J S Marron","Ming-Hui Chen","Yanyuan Ma","Cheng Yong Tang","Bing-Yi Jing","Michael Sherman","Rasmus Waagepetersen","Cun-Hui Zhang","Hongzhe Li","Jian Huang","David L Donoho","Grace Wahba","Annie Qu","Lawrence D Brown","Rui Song","Xiaotong Shen","Yichao Wu","Guang Cheng","Hui Zou","Trevor J Hastie","Zhezhen Jin","Howell Tong","Wang Zhou","Yong Zhou","Chenlei Leng","Yingcun Xia","Qingxia Chen","Byeong U Park","Abel Rodriguez","Bo Li","Kung-Sik Chan","Lan Zhang","Jacob Bien","Chiung-Yu Huang","Jiancheng Jiang","Daniela M Witten","Hansheng Wang","Yang Feng","Yacine Ait-Sahalia","Noelle I Samia","Dan Yu Lin","Hsin-Cheng Huang","Wei Pan","Hongtu Zhu","Joseph G Ibrahim","Donglin Zeng","Heng Peng","Runze Li","Yingying Fan","Jianqing Fan","Hua Liang","Tony Cai","Raymond J Carroll","Peter Hall","Hans-Georg Muller","Jing Qin"];
toshow_id=find(contains(author,key_author));
toshowg1=intersect(g1,toshow_id);
toshowg2=intersect(g2,toshow_id);
toshowg3=intersect(g3,toshow_id);
pU=cmdscale(squareform(pdist(tmp)),2);
figure;
plot(pU(g1,1),pU(g1,2),'.','color','b','MarkerSize',14)
hold on;
plot(pU(g2,1),pU(g2,2),'.','color','r','MarkerSize',14)
plot(pU(g3,1),pU(g3,2),'.','color','m','MarkerSize',14)
text(pU(toshowg1,1),pU(toshowg1,2),author(toshowg1),'color','b','fontsize',14)
text(pU(toshowg2,1),pU(toshowg2,2),author(toshowg2),'color','r','fontsize',14)
text(pU(toshowg3,1),pU(toshowg3,2),author(toshowg3),'color','m','fontsize',14)
set(gcf,'Position',[1000 593 900 750])
xlabel("PC 1 by SCORE and Multidimensional Scaling",'fontsize',14);
ylabel("PC 2 by SCORE and Multidimensional Scaling",'fontsize',14);

%%%%% SCORE and MDS and plot of authors relation -- alpha=0.0005
% tmpU=U{1};
% tmp=tmpU(:,2:k)./tmpU(:,1);tmp(isnan(tmp))=100;delta=15;tmp=trim(tmp,delta);
% [cidx, ctrs] = kmeans(tmp, 3, 'Replicates',20);
% g1=find(cidx==1);
% g2=find(cidx==2);
% g3=find(cidx==3);
% key_author=["David Dunson","Amy H Herring","Anirban Bhattacharya","Lawrence Carin","Yongtao Guan","Brian J Reich","Song Xi Chen","Fang Yao","Ying Wei","Heping Zhang","Harrison H Zhou","Debashis Paul","Yi Lin","Michael J Todd","Kani Chen","Qiwei Yao","Qi-Man Shao","Ming Yuan","Richard Samworth","J S Marron","Ming-Hui Chen","Yanyuan Ma","Cheng Yong Tang","Bing-Yi Jing","Michael Sherman","Rasmus Waagepetersen","Cun-Hui Zhang","Hongzhe Li","Jian Huang","David L Donoho","Grace Wahba","Annie Qu","Lawrence D Brown","Rui Song","Xiaotong Shen","Yichao Wu","Guang Cheng","Hui Zou","Trevor J Hastie","Zhezhen Jin","Howell Tong","Wang Zhou","Yong Zhou","Chenlei Leng","Yingcun Xia","Qingxia Chen","Byeong U Park","Abel Rodriguez","Bo Li","Kung-Sik Chan","Lan Zhang","Jacob Bien","Chiung-Yu Huang","Jiancheng Jiang","Daniela M Witten","Hansheng Wang","Yang Feng","Yacine Ait-Sahalia","Noelle I Samia","Dan Yu Lin","Hsin-Cheng Huang","Wei Pan","Hongtu Zhu","Joseph G Ibrahim","Donglin Zeng","Heng Peng","Runze Li","Yingying Fan","Jianqing Fan","Hua Liang","Tony Cai","Raymond J Carroll","Peter Hall","Hans-Georg Muller","Jing Qin"];
% toshow_id=find(contains(author,key_author));
% toshowg1=intersect(g1,toshow_id);
% toshowg2=intersect(g2,toshow_id);
% toshowg3=intersect(g3,toshow_id);
% pU=cmdscale(squareform(pdist(tmp)),2);
% figure;
% plot(pU(g1,1),pU(g1,2),'.','color','b','MarkerSize',14)
% hold on;
% plot(pU(g2,1),pU(g2,2),'.','color','r','MarkerSize',14)
% plot(pU(g3,1),pU(g3,2),'.','color','m','MarkerSize',14)
% text(pU(toshowg1,1),pU(toshowg1,2),author(toshowg1),'color','b','fontsize',14)
% text(pU(toshowg2,1),pU(toshowg2,2),author(toshowg2),'color','r','fontsize',14)
% text(pU(toshowg3,1),pU(toshowg3,2),author(toshowg3),'color','m','fontsize',14)
% set(gcf,'Position',[1000 593 900 750])
% xlabel("PC 1 by SCORE",'fontsize',14);
% ylabel("PC 2 by SCORE",'fontsize',14);





