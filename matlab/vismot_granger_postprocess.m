clear all;

cd /project/3011085.03/analysis/granger
load vismot_parcels

d=dir('*granger_pre.mat');
dr=dir('*granger_pre_timereversed.mat');
for k = 1:numel(d)
  load(d(k).name);
  G(:,:,:,k) = granger.grangerspctrm;
  load(dr(k).name);
  Gr(:,:,:,k) = granger.grangerspctrm;
end

% map onto 94 parcels
oldG = G;
oldGr = Gr;
clear G Gr;
for k = 1:size(oldG,4)
  for m = 1:size(oldG,3)
    tmp = oldG(:,:,m,k);
    tmp = P*tmp*P';
    G(:,:,m,k) = tmp;
    tmp = oldGr(:,:,m,k);
    tmp = P*tmp*P';
    Gr(:,:,m,k) = tmp;
  end
end
for k = 1:94
  G(k,k,:,:) = nan;
  Gr(k,k,:,:) = nan;
end

n = size(G,4);
design = [ones(1,n) ones(1,n)*2;1:n 1:n];

G(~isfinite(G))=0;
Gr(~isfinite(Gr)) = 0;

cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.computeprob = 'yes';
cfg.tail = 1;
stat = ft_statfun_depsamplesT(cfg,[reshape(G,[],19) reshape(Gr,[],19)],design);


%curr_dir = pwd;
%cd('~/matlab/fieldtrip/private');
%cd(curr_dir);
n = 94;
p = reshape(stat.prob, [n n 240]);

sP5  = sum(p<0.05,3)./240;
sP1  = sum(p<0.01,3)./240;
sP01 = sum(p<0.001,3)./240;
sP001 = sum(p<0.0001,3)./240;

ix1 = 1:47;
ix2 = 94:-1:48;

P5 = (sP5(ix1,ix1)+sP5(ix1,ix2)+sP5(ix2,ix1)+sP5(ix2,ix2))./4;
P1 = (sP1(ix1,ix1)+sP1(ix1,ix2)+sP1(ix2,ix1)+sP1(ix2,ix2))./4;
P01 = (sP01(ix1,ix1)+sP01(ix1,ix2)+sP01(ix2,ix1)+sP01(ix2,ix2))./4;
P001 = (sP001(ix1,ix1)+sP001(ix1,ix2)+sP001(ix2,ix1)+sP001(ix2,ix2))./4;

save('pthresholds_reversetime', 'sP5', 'sP1', 'sP01', 'sP001', 'P5', 'P1', 'P01', 'P001');

d=dir('*granger_pre.mat');
d1=dir('*prev1*');
d2=dir('*prev2*');
d3=dir('*prev3*');
d4=dir('*prev4*');
clear G;
for k = 1:numel(d1)
  k
  load(d(k).name);
  G(:,:,:,k) = granger.grangerspctrm;
  load(d1(k).name);
  G1(:,:,:,k) = granger.grangerspctrm;
  load(d2(k).name);
  G2(:,:,:,k) = granger.grangerspctrm;
  load(d3(k).name);
  G3(:,:,:,k) = granger.grangerspctrm;
  load(d4(k).name);
  G4(:,:,:,k) = granger.grangerspctrm;
end
for k = 1:size(G,4)
  for m = 1:size(G,3)
    tmp = G(:,:,m,k);
    tmp = P*tmp*P';
    newG(:,:,m,k) = tmp;
    tmp = G1(:,:,m,k);
    tmp = P*tmp*P';
    newG1(:,:,m,k) = tmp;
    tmp = G2(:,:,m,k);
    tmp = P*tmp*P';
    newG2(:,:,m,k) = tmp;
    tmp = G3(:,:,m,k);
    tmp = P*tmp*P';
    newG3(:,:,m,k) = tmp;
    tmp = G4(:,:,m,k);
    tmp = P*tmp*P';
    newG4(:,:,m,k) = tmp;
  end
end
for k = 1:94
  newG(k,k,:,:) = nan;
  newG1(k,k,:,:) = nan;
  newG2(k,k,:,:) = nan;
  newG3(k,k,:,:) = nan;
  newG4(k,k,:,:) = nan;
end
G = newG;
G1 = newG1; clear newG1;
G2 = newG2; clear newG2;
G3 = newG3; clear newG3;
G4 = newG4; clear newG4;

ix   = 47:-1:1;
mask = [P01>0 P01(:,ix)>0;P01(ix,:)>0 P01(ix,ix)>0];
%mask = [P001>0 P001>0;P001>0 P001>0];

dat = permute(G,[3 4 1 2]);
dat = reshape(dat,[240*19 n*n]);

option = [];
option.lambda = 0.2;
option.kernel = 'linear';

% ensure the mask is false on the diagonal
for k = 1:size(mask,1)
  mask(k,k) = false;
end
[A,S] = sparsenmfnnls_wrapper(10,option,dat(:,mask(:)));




for k = 1:240
  for m = 1:19
    tmp = G1(:,:,k,m);
    tmp(~mask) = nan;
    G1(:,:,k,m) = tmp;
    tmp = G2(:,:,k,m);
    tmp(~mask) = nan;
    G2(:,:,k,m) = tmp;
    tmp = G3(:,:,k,m);
    tmp(~mask) = nan;
    G3(:,:,k,m) = tmp;
    tmp = G4(:,:,k,m);
    tmp(~mask) = nan;
    G4(:,:,k,m) = tmp;
    tmp = G(:,:,k,m);
    tmp(~mask) = nan;
    G(:,:,k,m) = tmp;
    
  end
end

load vismot_parcels

for k = 1:numel(ulabel2)
  for m = 1:numel(ulabel2)
    selk = contains(ulabel,ulabel2{k});
    selm = contains(ulabel,ulabel2{m});
    G1x(k,m,:,:) = nanmedian(nanmedian(G1(selk,selm,:,:),1),2);
    G2x(k,m,:,:) = nanmedian(nanmedian(G2(selk,selm,:,:),1),2);
    G3x(k,m,:,:) = nanmedian(nanmedian(G3(selk,selm,:,:),1),2);
    G4x(k,m,:,:) = nanmedian(nanmedian(G4(selk,selm,:,:),1),2);
    Gx(k,m,:,:)  = nanmedian(nanmedian(G(selk,selm,:,:),1),2);
  end
end
for k = 1:size(G1x,1)
  G1x(k,k,:,:) = nan;
  G2x(k,k,:,:) = nan;
  G3x(k,k,:,:) = nan;
  G4x(k,k,:,:) = nan;
  Gx(k,k,:,:)  = nan;
end

n = 19;
design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgx.ivar=1;
cfgx.uvar=2;
stat13=ft_statfun_wilcoxon(cfgx,[reshape(G1x,[],19) reshape(G3x,[],19)],design);
stat42=ft_statfun_wilcoxon(cfgx,[reshape(G4x,[],19) reshape(G2x,[],19)],design);

T13 = reshape(stat13.stat,[16 16 240]);
T42 = reshape(stat42.stat,[16 16 240]);

re_indx = 16:-1:1;
T       = (T13+T42(re_indx,re_indx,:))./sqrt(2);

% further averaging into 6x6
pmat = blkdiag([1 1 1],[1 1],[1 1 1],[1 1 1],[1 1],[1 1 1]);
for k = 1:size(pmat,1)
  for m = 1:size(pmat,1)
    selk = pmat(k,:)>0;
    selm = pmat(m,:)>0;
    G1y(k,m,:,:) = nanmean(nanmean(G1x(selk,selm,:,:)));
    G2y(k,m,:,:) = nanmean(nanmean(G2x(selk,selm,:,:)));
    G3y(k,m,:,:) = nanmean(nanmean(G3x(selk,selm,:,:)));
    G4y(k,m,:,:) = nanmean(nanmean(G4x(selk,selm,:,:)));
  end
end
for k = 1:size(G1y,1)
  G1y(k,k,:,:) = nan;
  G2y(k,k,:,:) = nan;
  G3y(k,k,:,:) = nan;
  G4y(k,k,:,:) = nan;
end
re_indx = 6:-1:1;


cfgx.numrandomization = 1000;
cfgx.statistic = 'ft_statfun_wilcoxon';
stat13y = ft_statistics_montecarlo(cfgx,[reshape(G1y,[],19) reshape(G3y,[],19)],design);
stat42y = ft_statistics_montecarlo(cfgx,[reshape(G4y,[],19) reshape(G2y,[],19)],design);
staty   = ft_statistics_montecarlo(cfgx,[reshape(G1y+G4y(re_indx,re_indx,:,:),[],19) reshape(G3y+G2y(re_indx,re_indx,:,:),[],19)],design);

findx = nearest(granger.freq,8):nearest(granger.freq,12);
dat   = [reshape(mean(G1y(:,:,findx,:)+G4y(re_indx,re_indx,findx,:),3),[],19) reshape(mean(G3y(:,:,findx,:)+G2y(re_indx,re_indx,findx,:),3),[],19)];
stat  = ft_statistics_montecarlo(cfgx, dat, design);

