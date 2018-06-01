function [sx,sx1,sx2,allinside,dim,alpha1,alpha2,beta1,beta2] = doSourcestatisticsFastSlowCoh 

subjinfo;
names = {SUBJ(:).name};
names = names([1:3 5 6 8:11 13:16 18:20]);

dat1=zeros(16,33480,10,8);                                    
dat2=zeros(16,33480,10,8);                                    
dum=zeros(33480,1);
allinside=zeros(33480,1);
for k = 1:numel(names)
  load([names{k},'cohFastSlowSmooth'],'s');
  tmpin = s{1}.inside;
  for kk = 1:4
    tmp   = s{kk};
    num   = atanh(abs(tmp.coh1))-1./(tmp.df1-2);
    spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh1 = num;
    num   = atanh(abs(tmp.coh2))-1./(tmp.df2-2);
    spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh2 = num;
    num   = atanh(abs(tmp.coh3))-1./(tmp.df3-2);
    spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh3 = num;
    num   = atanh(abs(tmp.coh4))-1./(tmp.df4-2);
    spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh4 = num;
    s{kk} = tmp;
  end
  dat1(k,:,:,1)=s{1}.coh1;
  dat1(k,:,:,2)=s{2}.coh1;
  dat1(k,:,:,3)=s{3}.coh1;
  dat1(k,:,:,4)=s{4}.coh1;
  dat1(k,:,:,5)=s{1}.coh2;
  dat1(k,:,:,6)=s{2}.coh2;
  dat1(k,:,:,7)=s{3}.coh2;
  dat1(k,:,:,8)=s{4}.coh2;
  dat2(k,:,:,1)=s{1}.coh3;
  dat2(k,:,:,2)=s{2}.coh3;
  dat2(k,:,:,3)=s{3}.coh3;
  dat2(k,:,:,4)=s{4}.coh3;
  dat2(k,:,:,5)=s{1}.coh4;
  dat2(k,:,:,6)=s{2}.coh4;
  dat2(k,:,:,7)=s{3}.coh4;
  dat2(k,:,:,8)=s{4}.coh4;
  rt(k,1:2) = s{1}.rt;
  rt(k,3:4) = s{2}.rt;
  rt(k,5:6) = s{3}.rt;
  rt(k,7:8) = s{4}.rt;
  dum(s{1}.roi1) = dum(s{1}.roi1)+1;
  dum(s{1}.roi2) = dum(s{2}.roi2)-1;
  allinside(s{1}.inside) = allinside(tmpin)+1;
end
allinside=find(allinside==16);
ninside=numel(allinside);
dat1=permute(dat1(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
dat1=reshape(dat1,[ninside*10 16*8]);
dat2=permute(dat2(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
dat2=reshape(dat2,[ninside*10 16*8]);
design(1,:) = reshape(repmat(1:16, [8 1]), [1 128]);
design(2,:) = repmat([1 1 1 1 2 2 2 2], [1 16]);
design(3,:) = repmat([1 1 2 2 1 1 2 2], [1 16]);
design(4,:) = repmat([1 2 2 1 1 2 2 1], [1 16]);
cfg.ivar = [2 3 4];
cfg.uvar = 1;

dat1 = reshape(dat1, [ninside 10 16*8]);
dat2 = reshape(dat2, [ninside 10 16*8]);

alphadat1 = squeeze(mean(dat1(:,1:3,:),2));
alphadat2 = squeeze(mean(dat2(:,1:3,:),2));
betadat1 = squeeze(mean(dat1(:,5:8,:),2));
betadat2 = squeeze(mean(dat2(:,5:8,:),2));

dim    = pos2dim3d(s{1}.pos);
dummy  = zeros(dim)-2;
alpha1 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
alpha2= struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
beta1 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
beta2 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);

cfg.f = 'omnibus';
alpha1.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'main1';
alpha1.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'main2';
alpha1.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'main3';
alpha1.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'interaction12';
alpha1.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'interaction13';
alpha1.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'interaction23';
alpha1.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
cfg.f = 'interaction123';
alpha1.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');

M = blkdiag(ones(1,4),ones(1,4));
M = blkdiag(M,M);
M = blkdiag(M,M);
M = blkdiag(M,M);
M = blkdiag(M,M)./4;

malpha1 = ((alphadat1))*M';
malpha2 = ((alphadat2))*M';
mbeta1  = ((betadat1))*M';
mbeta2  = ((betadat2))*M';

% volumetric correction term
ca1 = 0.5.*mean(malpha1(:,1:16))-mean(malpha1(:,17:32));
ca2 = 0.5.*mean(malpha2(:,1:16))-mean(malpha2(:,17:32));
malpha1(:,1:16) = malpha1(:,1:16) - repmat(ca1,[size(malpha1,1) 1]);
malpha2(:,1:16) = malpha2(:,1:16) - repmat(ca2,[size(malpha2,1) 1]);
malpha1(:,17:32) = malpha1(:,17:32) - repmat(ca1,[size(malpha1,1) 1]);
malpha2(:,17:32) = malpha2(:,17:32) - repmat(ca2,[size(malpha2,1) 1]);

inside = zeros(dim);
inside(allinside) = 1;
inside = flipdim(inside,1);
inside(allinside) = inside(allinside)+1;
newinside = find(inside==2);
newoutside = find(inside<2);

[int,ia,ib] = intersect(allinside, newinside);
malpha1 = malpha1(ia,:);
malpha2 = malpha2(ia,:);
mbeta1 = mbeta1(ia,:);
mbeta2 = mbeta2(ia,:);
allinside = newinside;
alloutside = newoutside;

cfg2 = [];
cfg2.ivar = 1;
cfg2.uvar = 2;
cfg2.statistic = 'depsamplesT';
cfg2.numrandomization = 1000;
cfg2.correctm = 'fdr';
design2 = repmat([1 2],[1 16]);
design2(2,1:2:end) = 1:16;
design2(2,2:2:end) = 1:16;
sx1 = statistics_montecarlo(cfg2, malpha1, design2);
sx2 = statistics_montecarlo(cfg2, malpha2, design2);

%flip and pool
palpha = malpha1;
pbeta  = mbeta1;
dumvar = zeros(dim);
outvol = zeros(dim);
outvol(allinside) = outvol(allinside)+1;
outvol = flipdim(outvol,1);
outvol(allinside) = outvol(allinside)+1;
alloutside = find(outvol<2);

for k = 1:32
  dumvar(:) = 0;
  dumvar(allinside) = malpha2(:,k);
  dumvar = flipdim(dumvar,1);
  dumvar(newoutside) = 0;
  palpha(:,k) = (palpha(:,k)+dumvar(allinside))./sqrt(2);
  dumvar(:) = 0;
  dumvar(allinside) = mbeta2(:,k);
  dumvar = flipdim(dumvar,1);
  dumvar(newoutside) = 0;
  pbeta(:,k) = (pbeta(:,k)+dumvar(allinside))./sqrt(2);
  dumvar(:) = 0;
end
cfg2.numrandomization = 5000;
sx = statistics_montecarlo(cfg2, real(palpha), design2);
