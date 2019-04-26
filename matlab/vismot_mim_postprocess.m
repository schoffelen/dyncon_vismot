clear all;

load vismot_parcels
cd /project/3011085.03/analysis/mim


% if further dimension reduction to 6x6 (instead of (16x16)
if ~exist('dim','var'), dim=16; end % or 6
if dim==6
  load sixteen2six
  P3 = sixteen2six;
else
  P3=1;
end
freqs = 0.5:0.5:120;
nfreq = numel(freqs);
ncond=5;
n=19;

d = dir('*mim_pre.mat');
Mtmp2 = zeros(dim,dim,nfreq,ncond);
M_alltmp2 = zeros(dim,dim,nfreq);
for k = 1:numel(d)
  k
  load(d(k).name);
  Mtmp = cat(4,mim.mimspctrm);
  for m = 1:(nfreq*ncond)
    tmp = Mtmp(:,:,m);
    Mtmp2(:,:,m) = P3*P2*P*tmp*P'*P2'*P3';
  end
  M(:,:,:,:,k) = Mtmp2;
  
  M_alltmp = mim_all.mimspctrm;
  for m=1:nfreq
    tmp_all = M_alltmp(:,:,m);
    M_alltmp2(:,:,m) = P3*P2*P*tmp_all*P'*P2'*P3';
  end
  M_all(:,:,:,k) = M_alltmp2;
  
end
M = reshape(M,[dim dim nfreq ncond n]);


n = 19;
design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgx.ivar=1;
cfgx.uvar=2;
stat13=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,1,:),[],19) reshape(M(:,:,:,3,:),[],19)],design);
stat42=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,4,:),[],19) reshape(M(:,:,:,2,:),[],19)],design);
T13 = reshape(stat13.stat,[16 16 240]);
T42 = reshape(stat42.stat,[16 16 240]);

re_indx = 16:-1:1;
dat     = [reshape(M(:,:,:,1,:)+M(re_indx,re_indx,:,4,:),[],19) reshape(M(:,:,:,3,:)+M(re_indx,re_indx,:,2,:),[],19)];
stat    = ft_statfun_wilcoxon(cfgx, dat, design);
T       = reshape(stat.stat,[16 16 240]);


T13 = reshape(stat13.stat,[16 16 240]);
T42 = reshape(stat42.stat,[16 16 240]);
T   = reshape(stat.stat,  [16 16 240]);

stat15=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,1,:),[],19) reshape(M(:,:,:,5,:),[],19)],design);
stat25=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,2,:),[],19) reshape(M(:,:,:,5,:),[],19)],design);
stat35=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,3,:),[],19) reshape(M(:,:,:,5,:),[],19)],design);
stat45=ft_statfun_wilcoxon(cfgx,[reshape(M(:,:,:,4,:),[],19) reshape(M(:,:,:,5,:),[],19)],design);
T15   = reshape(stat15.stat,  [16 16 240]);
T25   = reshape(stat25.stat,  [16 16 240]);
T35   = reshape(stat35.stat,  [16 16 240]);
T45   = reshape(stat45.stat,  [16 16 240]);

cfgx.numrandomization = 1000;
cfgx.statistic = 'ft_statfun_wilcoxon';

foi   = 0:0.5:119.5;
stat  = ft_statistics_montecarlo(cfgx, dat, design);
T     = reshape(stat.stat,  [16 16 240]);

findx   = nearest(foi,8):nearest(foi,14);
alpha   = [reshape(mean(M(:,:,findx,1,:)+M(re_indx,re_indx,findx,4,:),3),[],19) reshape(mean(M(:,:,findx,3,:)+M(re_indx,re_indx,findx,2,:),3),[],19)];
statalpha = ft_statistics_montecarlo(cfgx, alpha, design);
findx   = nearest(foi,15):nearest(foi,30);
beta   = [reshape(mean(M(:,:,findx,1,:)+M(re_indx,re_indx,findx,4,:),3),[],19) reshape(mean(M(:,:,findx,3,:)+M(re_indx,re_indx,findx,2,:),3),[],19)];
statbeta = ft_statistics_montecarlo(cfgx, beta, design);
findx   = nearest(foi,30):nearest(foi,50);
gamma1   = [reshape(mean(M(:,:,findx,1,:)+M(re_indx,re_indx,findx,4,:),3),[],19) reshape(mean(M(:,:,findx,3,:)+M(re_indx,re_indx,findx,2,:),3),[],19)];
statgamma1 = ft_statistics_montecarlo(cfgx, gamma1, design);
findx   = nearest(foi,50):nearest(foi,70);
gamma2   = [reshape(mean(M(:,:,findx,1,:)+M(re_indx,re_indx,findx,4,:),3),[],19) reshape(mean(M(:,:,findx,3,:)+M(re_indx,re_indx,findx,2,:),3),[],19)];
statgamma2 = ft_statistics_montecarlo(cfgx, gamma2, design);
save('groupresults_freqbands', 'statalpha', 'statbeta', 'statgamma1', 'statgamma2');

