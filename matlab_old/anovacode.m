function [sx,sx1,sx2,allinside,dim,x,y,alpha1,alpha2,beta1,beta2] = anovacode

%dat=zeros(16,33480,10,8);                                    
%dum=zeros(33480,1);
%for k = 1:numel(names)
%  load(names{k},'s');
%  %dat(k,:,:,1)=sqrt(-(s{1}.df3-2)*log(1-double(s{1}.coh3).^2));
%  %dat(k,:,:,2)=sqrt(-(s{1}.df4-2)*log(1-double(s{1}.coh4).^2));
%  %dat(k,:,:,3)=sqrt(-(s{2}.df3-2)*log(1-double(s{2}.coh3).^2));
%  %dat(k,:,:,4)=sqrt(-(s{2}.df4-2)*log(1-double(s{2}.coh4).^2));
%  %dat(k,:,:,5)=sqrt(-(s{3}.df3-2)*log(1-double(s{3}.coh3).^2));
%  %dat(k,:,:,6)=sqrt(-(s{3}.df4-2)*log(1-double(s{3}.coh4).^2));
%  %dat(k,:,:,7)=sqrt(-(s{4}.df3-2)*log(1-double(s{4}.coh3).^2));
%  %dat(k,:,:,8)=sqrt(-(s{4}.df4-2)*log(1-double(s{4}.coh4).^2));
%  dat(k,:,:,1)=sqrt(-(s{1}.df1-2)*log(1-double(s{1}.coh1).^2));
%  dat(k,:,:,2)=sqrt(-(s{1}.df2-2)*log(1-double(s{1}.coh2).^2));
%  dat(k,:,:,3)=sqrt(-(s{2}.df1-2)*log(1-double(s{2}.coh1).^2));
%  dat(k,:,:,4)=sqrt(-(s{2}.df2-2)*log(1-double(s{2}.coh2).^2));
%  dat(k,:,:,5)=sqrt(-(s{3}.df1-2)*log(1-double(s{3}.coh1).^2));
%  dat(k,:,:,6)=sqrt(-(s{3}.df2-2)*log(1-double(s{3}.coh2).^2));
%  dat(k,:,:,7)=sqrt(-(s{4}.df1-2)*log(1-double(s{4}.coh1).^2));
%  dat(k,:,:,8)=sqrt(-(s{4}.df2-2)*log(1-double(s{4}.coh2).^2));
%  rt(k,1:2) = s{1}.rt;
%  rt(k,3:4) = s{2}.rt;
%  rt(k,5:6) = s{3}.rt;
%  rt(k,7:8) = s{4}.rt;
%  dum(s{1}.roi1) = dum(s{1}.roi1)+1;
%  dum(s{1}.roi2) = dum(s{2}.roi2)-1;
%end
%dat=permute(dat,[2 3 4 1]); %put subjects at slowest changing dimension
%dat=reshape(dat,[33480*10 16*8]);
%
%cfg.f = 'main1';
%stat1 = statfun_anova2x2rm(cfg, dat, design);
%cfg.f = 'main2';
%stat2 = statfun_anova2x2rm(cfg, dat, design);
%cfg.f = 'interaction';
%stat3 = statfun_anova2x2rm(cfg, dat, design);
%
%stat1.stat = reshape(stat1.stat,[33480 10]);     
%stat2.stat = reshape(stat2.stat,[33480 10]);     
%stat3.stat = reshape(stat3.stat,[33480 10]);     
% 
%dat=zeros(33480,10,16);
%alldat=zeros(33480,10,8,16);
%dum=zeros(33480,1);
%alldof=[];
%for k = 1:numel(names)
%  tmp = zeros(33480,10,8);
%  load(names{k},'s');
%  %tmp(:,:,1)=double(s{1}.coh3);
%  %tmp(:,:,2)=double(s{1}.coh4);
%  %tmp(:,:,3)=double(s{2}.coh3);
%  %tmp(:,:,4)=double(s{2}.coh4);
%  %tmp(:,:,5)=double(s{3}.coh3);
%  %tmp(:,:,6)=double(s{3}.coh4);
%  %tmp(:,:,7)=double(s{4}.coh3);
%  %tmp(:,:,8)=double(s{4}.coh4);
%  tmp(:,:,1)=double(s{1}.coh1);
%  tmp(:,:,2)=double(s{1}.coh2);
%  tmp(:,:,3)=double(s{2}.coh1);
%  tmp(:,:,4)=double(s{2}.coh2);
%  tmp(:,:,5)=double(s{3}.coh1);
%  tmp(:,:,6)=double(s{3}.coh2);
%  tmp(:,:,7)=double(s{4}.coh1);
%  tmp(:,:,8)=double(s{4}.coh2);
%  dof(1)=s{1}.df1;
%  dof(2)=s{1}.df2;
%  dof(3)=s{2}.df1;
%  dof(4)=s{2}.df2;
%  dof(5)=s{3}.df1;
%  dof(6)=s{3}.df2;
%  dof(7)=s{4}.df1;
%  dof(8)=s{4}.df2;
%  chi=coh2chi2(tmp,dof);
%  dat(:,:,k)=chi;
%  alldof=[alldof dof];
%  alldat(:,:,:,k)=tmp;
%  clear chi dof tmp;
%end
%
%dat1=zeros(16,33480,10,4);                                    
%dat2=zeros(16,33480,10,4);                                    
%dum=zeros(33480,1);
%allinside=zeros(33480,1);
%for k = 1:numel(names)
%  load(names{k},'s');
%  tmpin = s{1}.inside;
%  for kk = 1:4
%    tmp   = s{kk};
%    num   = atanh(abs(tmp.coh1))-atanh(abs(tmp.coh2));
%    denom = sqrt(1./(tmp.df1-2) + 1./(tmp.df2-2));
%    dummy = num./denom;
%    spm_smooth(dummy,dummy,1.5);dummy(tmpin)=dummy(tmpin)-nanmean(dummy(tmpin));
%    tmp.dcoh1 = dummy;
%    num   = atanh(abs(tmp.coh3))-atanh(abs(tmp.coh4));
%    denom = sqrt(1./(tmp.df3-2) + 1./(tmp.df4-2));
%    dummy = num./denom;
%    spm_smooth(dummy,dummy,1.5);dummy(tmpin)=dummy(tmpin)-nanmean(dummy(tmpin));
%    tmp.dcoh2 = dummy;
%    s{kk} = tmp;
%  end
%  dat1(k,:,:,1)=s{1}.dcoh1;
%  dat1(k,:,:,2)=s{2}.dcoh1;
%  dat1(k,:,:,3)=s{3}.dcoh1;
%  dat1(k,:,:,4)=s{4}.dcoh1;
%  dat2(k,:,:,1)=s{1}.dcoh2;
%  dat2(k,:,:,2)=s{2}.dcoh2;
%  dat2(k,:,:,3)=s{3}.dcoh2;
%  dat2(k,:,:,4)=s{4}.dcoh2;
%  rt(k,1:2) = s{1}.rt;
%  rt(k,3:4) = s{2}.rt;
%  rt(k,5:6) = s{3}.rt;
%  rt(k,7:8) = s{4}.rt;
%  dum(s{1}.roi1) = dum(s{1}.roi1)+1;
%  dum(s{1}.roi2) = dum(s{2}.roi2)-1;
%  allinside(s{1}.inside) = allinside(tmpin)+1;
%end
%allinside=find(allinside==16);
%ninside=numel(allinside);
%dat1=permute(dat1(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
%dat1=reshape(dat1,[ninside*10 16*4]);
%dat2=permute(dat2(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
%dat2=reshape(dat2,[ninside*10 16*4]);
%design(1,:) = reshape(repmat(1:16, [4 1]), [1 64]);
%design(2,:) = repmat([1 1 2 2], [1 16]);
%design(3,:) = repmat([1 2 2 1], [1 16]);
%cfg.ivar = [2 3];
%cfg.uvar = 1;
%
%M=ones(1,4).*0.25;
%M=blkdiag(M,M);
%M=blkdiag(M,M);
%M=blkdiag(M,M);
%M=blkdiag(M,M)';
%mdat1 = dat1*M;
%mdat2 = dat2*M;
%dat1 = dat1 - mdat1(:,design(1,:));
%dat2 = dat2 - mdat2(:,design(1,:));
%
%cfg.f = 'main1';
%stat1 = statfun_anova2x2rm(cfg, dat1, design);
%cfg.f = 'main2';
%stat2 = statfun_anova2x2rm(cfg, dat1, design);
%cfg.f = 'interaction';
%stat3 = statfun_anova2x2rm(cfg, dat1, design);
%x = s{1};
%x.stat1 = zeros(33480,10);
%x.stat2 = zeros(33480,10);
%x.stat3 = zeros(33480,10);
%x.stat1(allinside,:) = reshape(stat1.stat, [ninside 10]);
%x.stat2(allinside,:) = reshape(stat2.stat, [ninside 10]);
%x.stat3(allinside,:) = reshape(stat3.stat, [ninside 10]);
%
%cfg.f = 'main1';
%stat1 = statfun_anova2x2rm(cfg, dat2, design);
%cfg.f = 'main2';
%stat2 = statfun_anova2x2rm(cfg, dat2, design);
%cfg.f = 'interaction';
%stat3 = statfun_anova2x2rm(cfg, dat2, design);
%y = s{1};
%y.stat1 = zeros(33480,10);
%y.stat2 = zeros(33480,10);
%y.stat3 = zeros(33480,10);
%y.stat1(allinside,:) = reshape(stat1.stat, [ninside 10]);
%y.stat2(allinside,:) = reshape(stat2.stat, [ninside 10]);
%y.stat3(allinside,:) = reshape(stat3.stat, [ninside 10]);

subjinfo;
names = {SUBJ(:).name};
names = names([1:3 5 6 8:11 13:16 18:20]);

dat1=zeros(16,33480,31,8);                                    
dat2=zeros(16,33480,31,8);                                    
dum=zeros(33480,1);
allinside=zeros(33480,1);
for k = 1:numel(names)
  load([names{k},'cohFastSlowSmooth'],'s');
  tmpin = s{1}.inside;
  for kk = 1:4
    tmp   = s{kk};
    num   = atanh(abs(tmp.coh1))-1./(tmp.df1-2);
    %spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh1 = num;
    num   = atanh(abs(tmp.coh2))-1./(tmp.df2-2);
    %spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh2 = num;
    num   = atanh(abs(tmp.coh3))-1./(tmp.df3-2);
    %spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh3 = num;
    num   = atanh(abs(tmp.coh4))-1./(tmp.df4-2);
    %spm_smooth(num,num,1.5);%num(tmpin)=num(tmpin)-nanmean(num(tmpin));
    tmp.coh4 = num;
    s{kk} = tmp;
  end
  dat1(k,:,:,1)=s{1}.coh1; %ROI left, RT fast,  condition Left  Congruent  
  dat1(k,:,:,2)=s{2}.coh1; %ROI left, RT fast,  condition Left  Incongruent
  dat1(k,:,:,3)=s{3}.coh1; %ROI left, RT fast,  condition Right Congruent
  dat1(k,:,:,4)=s{4}.coh1; %ROI left, RT fast,  condition Right Incongruent  
  dat1(k,:,:,5)=s{1}.coh2; %ROI left, RT slow,  condition Left  Congruent  
  dat1(k,:,:,6)=s{2}.coh2; %ROI left, RT slow,  condition Left  Incongruent
  dat1(k,:,:,7)=s{3}.coh2; %ROI left, RT slow,  condition Right Congruent
  dat1(k,:,:,8)=s{4}.coh2; %ROI left, RT slow,  condition Right Incongruent  
  dat2(k,:,:,1)=s{1}.coh3; %ROI right, RT fast, condition Left  Congruent  
  dat2(k,:,:,2)=s{2}.coh3; %ROI right, RT fast, condition Left  Incongruent
  dat2(k,:,:,3)=s{3}.coh3; %ROI right, RT fast, condition Right Congruent
  dat2(k,:,:,4)=s{4}.coh3; %ROI right, RT fast, condition Right Inongruent  
  dat2(k,:,:,5)=s{1}.coh4; %ROI right, RT slow, condition Left  Congruent  
  dat2(k,:,:,6)=s{2}.coh4; %ROI right, RT slow, condition Left  Incongruent
  dat2(k,:,:,7)=s{3}.coh4; %ROI right, RT slow, condition Right Congruent
  dat2(k,:,:,8)=s{4}.coh4; %ROI right, RT slow, condition Right Incongruent  
  rt(k,1:2) = s{1}.rt;
  rt(k,3:4) = s{2}.rt;
  rt(k,5:6) = s{3}.rt;
  rt(k,7:8) = s{4}.rt;
  dof(k,1)  = s{1}.df1;
  dof(k,2)  = s{2}.df1;
  dof(k,3)  = s{3}.df1;
  dof(k,4)  = s{4}.df1;
  dum(s{1}.roi1) = dum(s{1}.roi1)+1;
  dum(s{1}.roi2) = dum(s{2}.roi2)-1;
  allinside(tmpin) = allinside(tmpin)+1;
end
rt = rt(:,[1 3 5 7 2 4 6 8]);
allinside=find(allinside==16);
ninside=numel(allinside);
dat1=permute(dat1(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
dat1=reshape(dat1,[ninside*31 16*8]);
dat2=permute(dat2(:,allinside,:,:),[2 3 4 1]); %put subjects at slowest changing dimension
dat2=reshape(dat2,[ninside*31 16*8]);
design(1,:) = reshape(repmat(1:16, [8 1]), [1 128]); 
design(2,:) = repmat([1 1 1 1 2 2 2 2], [1 16]); %Fast vs Slow
design(3,:) = repmat([1 1 2 2 1 1 2 2], [1 16]); %Condition left vs condition right
design(4,:) = repmat([1 2 2 1 1 2 2 1], [1 16]); %Stim in right visual hemi vs stim in left visual hemi
cfg.ivar = [2 3 4];
cfg.uvar = 1;

%M=ones(1,8)./8;
%M=blkdiag(M,M);
%M=blkdiag(M,M);
%M=blkdiag(M,M);
%M=blkdiag(M,M)';
%mdat1 = dat1*M;
%mdat2 = dat2*M;
%dat1 = dat1 - mdat1(:,design(1,:));
%dat2 = dat2 - mdat2(:,design(1,:));

cfg.f = 'main1';
stat1 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'main2';
stat2 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'main3';
stat3 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'interaction12';
stat12 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'interaction13';
stat13 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'interaction23';
stat23 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'interaction123';
stat123 = statfun_anova2x2x2rm(cfg, dat1, design);
cfg.f = 'omnibus';
stat0 = statfun_anova2x2x2rm(cfg, dat1, design);

x = s{1};
x.stat0 = zeros(33480,31);
x.stat1 = zeros(33480,31);
x.stat2 = zeros(33480,31);
x.stat3 = zeros(33480,31);
x.stat12 = zeros(33480,31);
x.stat13 = zeros(33480,31);
x.stat23 = zeros(33480,31);
x.stat123 = zeros(33480,31);
x.stat0(allinside,:) = reshape(stat0.stat, [ninside 31]);
x.stat1(allinside,:) = reshape(stat1.stat, [ninside 31]);
x.stat2(allinside,:) = reshape(stat2.stat, [ninside 31]);
x.stat3(allinside,:) = reshape(stat3.stat, [ninside 31]);
x.stat12(allinside,:) = reshape(stat12.stat, [ninside 31]);
x.stat13(allinside,:) = reshape(stat13.stat, [ninside 31]);
x.stat23(allinside,:) = reshape(stat23.stat, [ninside 31]);
x.stat123(allinside,:) = reshape(stat123.stat, [ninside 31]);

cfg.f = 'main1';
stat1 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'main2';
stat2 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'main3';
stat3 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'interaction12';
stat12 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'interaction13';
stat13 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'interaction23';
stat23 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'interaction123';
stat123 = statfun_anova2x2x2rm(cfg, dat2, design);
cfg.f = 'omnibus';
stat0 = statfun_anova2x2x2rm(cfg, dat2, design);

y = s{1};
y.stat0 = zeros(33480,31);
y.stat1 = zeros(33480,31);
y.stat2 = zeros(33480,31);
y.stat3 = zeros(33480,31);
y.stat12 = zeros(33480,31);
y.stat13 = zeros(33480,31);
y.stat23 = zeros(33480,31);
y.stat123 = zeros(33480,31);
y.stat0(allinside,:) = reshape(stat0.stat, [ninside 31]);
y.stat1(allinside,:) = reshape(stat1.stat, [ninside 31]);
y.stat2(allinside,:) = reshape(stat2.stat, [ninside 31]);
y.stat3(allinside,:) = reshape(stat3.stat, [ninside 31]);
y.stat12(allinside,:) = reshape(stat12.stat, [ninside 31]);
y.stat13(allinside,:) = reshape(stat13.stat, [ninside 31]);
y.stat23(allinside,:) = reshape(stat23.stat, [ninside 31]);
y.stat123(allinside,:) = reshape(stat123.stat, [ninside 31]);

dat1 = reshape(dat1, [ninside 31 16*8]);
dat2 = reshape(dat2, [ninside 31 16*8]);

alphadat1 = squeeze(mean(dat1(:,1:3,:),2));
alphadat2 = squeeze(mean(dat2(:,1:3,:),2));
betadat1  = squeeze(mean(dat1(:,5:10,:),2));
betadat2  = squeeze(mean(dat2(:,5:10,:),2));
gammadat1 = squeeze(mean(dat1(:,16:21,:),2));
gammadat2 = squeeze(mean(dat2(:,16:21,:),2));
gammadat1b = squeeze(mean(dat1(:,23:28,:),2));
gammadat2b = squeeze(mean(dat2(:,23:28,:),2));

dim    = pos2dim3d(y.pos);
dummy  = zeros(dim)-2;
alpha1 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
alpha2 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
beta1  = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
beta2  = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
gamma1 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
gamma2 = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
gamma1b = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);
gamma2b = struct('m0',dummy,'m1',dummy,'m2',dummy,'m3',dummy,'i12',dummy,'i13',dummy,'i23',dummy,'i123',dummy);

cfg.f = 'omnibus';
alpha1.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.m0(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'main1';
alpha1.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.m1(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'main2';
alpha1.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.m2(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'main3';
alpha1.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.m3(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'interaction12';
alpha1.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.i12(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'interaction13';
alpha1.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.i13(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'interaction23';
alpha1.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.i23(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');
cfg.f = 'interaction123';
alpha1.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat1, design), 'stat');
alpha2.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, alphadat2, design), 'stat');
 beta1.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat1, design), 'stat');
 beta2.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg, betadat2, design), 'stat');
gamma1.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1, design), 'stat');
gamma2.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2, design), 'stat');
gamma1b.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat1b, design), 'stat');
gamma2b.i123(allinside) = getfield(statfun_anova2x2x2rm(cfg,gammadat2b, design), 'stat');

%M = blkdiag(ones(1,4),ones(1,4));
%M = blkdiag(M,M);
%M = blkdiag(M,M);
%M = blkdiag(M,M);
%M = blkdiag(M,M)./4;

%M1 = M;
%M2 = M;

tmp = dof;
tmp(:,[1 3 4]) = 0;
tmp(:,2) = sqrt(tmp(:,2)-2)./2;
tmp = [tmp -tmp];
M1 = tmp(1,:);
for k = 2:size(tmp,1)
  M1 = blkdiag(M1,tmp(k,:));
end
tmp = dof;
tmp(:,[1 2 3]) = 0;
tmp(:,[4]) = sqrt(tmp(:,[4])-2)./2;
tmp = [tmp -tmp];
M2 = tmp(1,:);
for k = 2:size(tmp,1)
  M2 = blkdiag(M2,tmp(k,:));
end


malpha1 = ((alphadat1))*M1';
malpha2 = ((alphadat2))*M2';
mbeta1  = ((betadat1))*M1';
mbeta2  = ((betadat2))*M2';
mgamma1 = ((gammadat1))*M1';
mgamma2 = ((gammadat2))*M2';
mgamma1b  = ((gammadat1b))*M1';
mgamma2b  = ((gammadat2b))*M2';

%% volumetric correction term
%ca1 = 0.5.*(mean(malpha1(:,1:16))-mean(malpha1(:,17:32)));
%ca2 = 0.5.*(mean(malpha2(:,1:16))-mean(malpha2(:,17:32)));
%malpha1(:,1:16)  = malpha1(:,1:16) - repmat(ca1,[size(malpha1,1) 1]);
%malpha2(:,1:16)  = malpha2(:,1:16) - repmat(ca2,[size(malpha2,1) 1]);
%malpha1(:,17:32) = malpha1(:,17:32) - repmat(ca1,[size(malpha1,1) 1]);
%malpha2(:,17:32) = malpha2(:,17:32) - repmat(ca2,[size(malpha2,1) 1]);
%% volumetric correction term
%ca1 = 0.5.*(mean(mbeta1(:,1:16))-mean(mbeta1(:,17:32)));
%ca2 = 0.5.*(mean(mbeta2(:,1:16))-mean(mbeta2(:,17:32)));
%mbeta1(:,1:16)  = mbeta1(:,1:16) - repmat(ca1,[size(mbeta1,1) 1]);
%mbeta2(:,1:16)  = mbeta2(:,1:16) - repmat(ca2,[size(mbeta2,1) 1]);
%mbeta1(:,17:32) = mbeta1(:,17:32) - repmat(ca1,[size(mbeta1,1) 1]);
%mbeta2(:,17:32) = mbeta2(:,17:32) - repmat(ca2,[size(mbeta2,1) 1]);
%% volumetric correction term
%ca1 = 0.5.*(mean(mgamma1(:,1:16))-mean(mgamma1(:,17:32)));
%ca2 = 0.5.*(mean(mgamma2(:,1:16))-mean(mgamma2(:,17:32)));
%mgamma1(:,1:16)  = mgamma1(:,1:16) - repmat(ca1,[size(mgamma1,1) 1]);
%mgamma2(:,1:16)  = mgamma2(:,1:16) - repmat(ca2,[size(mgamma2,1) 1]);
%mgamma1(:,17:32) = mgamma1(:,17:32) - repmat(ca1,[size(mgamma1,1) 1]);
%mgamma2(:,17:32) = mgamma2(:,17:32) - repmat(ca2,[size(mgamma2,1) 1]);
%% volumetric correction term
%ca1 = 0.5.*(mean(mgamma1b(:,1:16))-mean(mgamma1b(:,17:32)));
%ca2 = 0.5.*(mean(mgamma2b(:,1:16))-mean(mgamma2b(:,17:32)));
%mgamma1b(:,1:16)  = mgamma1b(:,1:16) - repmat(ca1,[size(mgamma1b,1) 1]);
%mgamma2b(:,1:16)  = mgamma2b(:,1:16) - repmat(ca2,[size(mgamma2b,1) 1]);
%mgamma1b(:,17:32) = mgamma1b(:,17:32) - repmat(ca1,[size(mgamma1b,1) 1]);
%mgamma2b(:,17:32) = mgamma2b(:,17:32) - repmat(ca2,[size(mgamma2b,1) 1]);
%%THE FOLLOWING DOES NOT MAKE SENSE BECAUSE WE ARE DEALING WITH COHERENCE
%%malpha1 = malpha1 - repmat(median(malpha1), [size(malpha1,1) 1]);
%%malpha2 = malpha2 - repmat(median(malpha2), [size(malpha2,1) 1]);
%%mbeta1 = mbeta1 - repmat(median(mbeta1), [size(mbeta1,1) 1]);
%%mbeta2 = mbeta2 - repmat(median(mbeta2), [size(mbeta2,1) 1]);
%%mgamma1 = mgamma1 - repmat(median(mgamma1), [size(mgamma1,1) 1]);
%%mgamma2 = mgamma2 - repmat(median(mgamma2), [size(mgamma2,1) 1]);
%%mgamma1b = mgamma1b - repmat(median(mgamma1b), [size(mgamma1b,1) 1]);
%%mgamma2b = mgamma2b - repmat(median(mgamma2b), [size(mgamma2b,1) 1]);
malpha1 = malpha1 - repmat(median(malpha1), [size(malpha1,1) 1]);
malpha2 = malpha2 - repmat(median(malpha2), [size(malpha2,1) 1]);
mbeta1 = mbeta1 - repmat(median(mbeta1), [size(mbeta1,1) 1]);
mbeta2 = mbeta2 - repmat(median(mbeta2), [size(mbeta2,1) 1]);
mgamma1 = mgamma1 - repmat(median(mgamma1), [size(mgamma1,1) 1]);
mgamma2 = mgamma2 - repmat(median(mgamma2), [size(mgamma2,1) 1]);
mgamma1b = mgamma1b - repmat(median(mgamma1b), [size(mgamma1b,1) 1]);
mgamma2b = mgamma2b - repmat(median(mgamma2b), [size(mgamma2b,1) 1]);

nsubj = size(malpha1,2);
malpha1(:,nsubj+(1:nsubj)) = 0;
malpha2(:,nsubj+(1:nsubj)) = 0;
mbeta1(:,nsubj+(1:nsubj)) = 0;
mbeta2(:,nsubj+(1:nsubj)) = 0;
mgamma1(:,nsubj+(1:nsubj)) = 0;
mgamma2(:,nsubj+(1:nsubj)) = 0;
mgamma1b(:,nsubj+(1:nsubj)) = 0;
mgamma2b(:,nsubj+(1:nsubj)) = 0;

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
mgamma1 = mgamma1(ia,:);
mgamma2 = mgamma2(ia,:);
mgamma1b = mgamma1b(ia,:);
mgamma2b = mgamma2b(ia,:);
allinside = newinside;
alloutside = newoutside;

cfg2 = [];
cfg2.ivar = 1;
cfg2.uvar = 2;
cfg2.statistic = 'depsamplesT';
cfg2.numrandomization = 10;
cfg2.correctm = 'fdr';
%design2 = repmat([1 2],[1 16]);
%design2(2,1:2:end) = 1:16;
%design2(2,2:2:end) = 1:16;
design2 = [ones(1,16) ones(1,16)*2; 1:16 1:16];
sax1 = statistics_montecarlo(cfg2, malpha1, design2);
sax2 = statistics_montecarlo(cfg2, malpha2, design2);
sbx1 = statistics_montecarlo(cfg2, mbeta1, design2);
sbx2 = statistics_montecarlo(cfg2, mbeta2, design2);
sgx1 = statistics_montecarlo(cfg2, mgamma1, design2);
sgx2 = statistics_montecarlo(cfg2, mgamma2, design2);
sgx1b = statistics_montecarlo(cfg2, mgamma1b, design2);
sgx2b = statistics_montecarlo(cfg2, mgamma2b, design2);
:w
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
  tmp = [];
  dumvar(:) = 0;
  dumvar(allinside) = malpha1(:,k);
  dumvar(alloutside) = 0;
  spm_smooth(dumvar, dumvar, 1.5);
  tmp = dumvar;
  dumvar(:) = 0;
  dumvar(allinside) = malpha2(:,k);
  dumvar = flipdim(dumvar,1);
  dumvar(alloutside) = 0;
  spm_smooth(dumvar, dumvar, 1.5);
  palpha(:,k) = (tmp(allinside)+dumvar(allinside))./sqrt(2);
  tmp = [];
  dumvar(:) = 0;
  dumvar(allinside) = mbeta1(:,k);
  dumvar(alloutside) = 0;
  spm_smooth(dumvar, dumvar, 1.5);
  tmp = dumvar;
  dumvar(:) = 0;
  dumvar(allinside) = mbeta2(:,k);
  dumvar = flipdim(dumvar,1);
  dumvar(alloutside) = 0;
  spm_smooth(dumvar, dumvar, 1.5);
  pbeta(:,k) = (tmp(allinside)+dumvar(allinside))./sqrt(2);
end
cfg2.numrandomization = 1000;
cfg2.clusteralpha = 0.025;
cfg2.correctm = 'cluster';
cfg2.clusterthreshold = 'nonparametric_individual';
sx = statistics_montecarlo(cfg2, real(palpha), design2, 'issource', 1);
