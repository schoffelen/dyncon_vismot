
load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
frequency = [10 22 38 42 58 62 78 82];
n=numel(list);

toi = 'pre';
roi_to = 'roi';
load(fullfile([alldir, 'analysis/source/roi.mat']));
nregions = 2*(size(roi,1)-1);

cnt=1;
for k=2:size(foi,1)
  tmpfreqs = foi{k,4};
  for subj=1:n
    sz = numel(tmpfreqs);
    tmpzx13 = [];
    tmpzx42 = [];
    tmpcoh1 = [];
    tmpcoh2 = [];
    tmpcoh3 = [];
    tmpcoh4 = [];
    for f = 1:sz
      
      filename = fullfile([datadir, sprintf('%s_coh6d4mm_%s_roi2%s_%03d', list{subj},toi,roi_to, tmpfreqs(f))]);
      dum = load(filename);
      
      tmpzx13(:,:,f) = dum.zx13;
      tmpzx42(:,:,f) = dum.zx42;
      tmpcoh1(:,:,f) = abs(dum.coh(1).coh);
      tmpcoh2(:,:,f) = abs(dum.coh(2).coh);
      tmpcoh3(:,:,f) = abs(dum.coh(3).coh);
      tmpcoh4(:,:,f) = abs(dum.coh(4).coh);
    end
    zx13{cnt}(subj,:,:) = nanmean(tmpzx13,3);
    zx42{cnt}(subj,:,:) = nanmean(tmpzx42,3);
    coh1{cnt}(subj,:,:) = nanmean(tmpcoh1,3);
    coh2{cnt}(subj,:,:) = nanmean(tmpcoh2,3);
    coh3{cnt}(subj,:,:) = nanmean(tmpcoh3,3);
    coh4{cnt}(subj,:,:) = nanmean(tmpcoh4,3);
  end
  cnt=cnt+1;
end

% % pool zx13 and zx42: flip zx42 and average
% for k=1:numel(zx42)
%   nloc = size(zx42{k},3)/2;
%   zx42

% calculate average within/between hemisphere coherence
% seperately for conditions 1-3 and 4-2
for k=1:numel(zx13)
  nloc = size(zx13{k},3)/2;
  x = ones(2*nloc);
  for subj=1:n
    x13tmp = squeeze(zx13{k}(subj,:,:));
    x13tmp(find(triu(x))) = nan;
    
    zx13_between(subj,k) = nanmean(nanmean(x13tmp(nloc+1:2*nloc, 1:nloc)));
    zx13_within(subj,k) = nanmean([nanmean(nanmean(x13tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x13tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    
    x42tmp = squeeze(zx42{k}(subj,:,:));
    x42tmp(find(triu(x))) = nan;
    zx42_between(subj,k) = nanmean(nanmean(x42tmp(nloc+1:2*nloc, 1:nloc)));
    zx42_within(subj,k) = nanmean([nanmean(nanmean(x42tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x42tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    x1tmp = squeeze(coh1{k}(subj,:,:));
    x1tmp(find(triu(x))) = nan;
    coh1_between(subj,k) = nanmean(nanmean(x1tmp(nloc+1:2*nloc, 1:nloc)));
    coh1_within(subj,k) = nanmean([nanmean(nanmean(x1tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x1tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    x2tmp = squeeze(coh2{k}(subj,:,:));
    x2tmp(find(triu(x))) = nan;
    coh2_between(subj,k) = nanmean(nanmean(x2tmp(nloc+1:2*nloc, 1:nloc)));
    coh2_within(subj,k) = nanmean([nanmean(nanmean(x2tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x2tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    x3tmp = squeeze(coh3{k}(subj,:,:));
    x3tmp(find(triu(x))) = nan;
    coh3_between(subj,k) = nanmean(nanmean(x3tmp(nloc+1:2*nloc, 1:nloc)));
    coh3_within(subj,k) = nanmean([nanmean(nanmean(x3tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x3tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    x4tmp = squeeze(coh4{k}(subj,:,:));
    x4tmp(find(triu(x))) = nan;
    coh4_between(subj,k) = nanmean(nanmean(x4tmp(nloc+1:2*nloc, 1:nloc)));
    coh4_within(subj,k) = nanmean([nanmean(nanmean(x4tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x4tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
  end
end

% pool conditions 13 and 42
zx_between = (zx13_between + zx42_between)./2;
zx_within = (zx13_within + zx42_within)./2;

% pool conditions 1+4 and 1+3
coh_between_C = (coh1_between + coh4_between)./2;
coh_between_IC = (coh3_between + coh2_between)./2;
coh_within_C = (coh1_within + coh4_within)./2;
coh_within_IC = (coh3_within + coh2_within)./2;

% average over frequencies
zx_between = nanmean(zx_between,2);
zx_within = nanmean(zx_within,2);
coh_between_C = nanmean(coh_between_C, 2);
coh_between_IC = nanmean(coh_between_IC, 2);
coh_within_C = nanmean(coh_within_C, 2);
coh_within_IC = nanmean(coh_within_IC, 2);


% we expect within hemisphere connectivity to be higher for conditions 1-3
% and 4-2, w.r.t. between hemispheres connectivity: i.e.
% zx_within>zxbetween
[w.H,w.P,w.CI,w.STATS] = ttest(zx_within, zx_between, 'tail', 'right');
within_vs_between_hemiC = w;

[w.H,w.P,w.CI,w.STATS] = ttest(coh_within_C, coh_within_IC, 'tail', 'right');
within_CvsIC = w;

[w.H,w.P,w.CI,w.STATS] = ttest(coh_between_C, coh_between_IC, 'tail', 'left');
between_CvsIC = w;


%% Is there any connectivity pattern significant individually (i.e. larger/smaller for C than IC)?
% combine 1-3 and 4-2
zx42_flipped = zx42;
coh2_flipped = coh2;
coh4_flipped = coh4;
for k=1:numel(zx42)
  nloc = size(zx13{k},3)/2;
  zx42_flipped{k}(:,1:nloc,1:nloc) = zx42{k}(:,nloc+1:2*nloc, nloc+1:2*nloc);
  zx42_flipped{k}(:,nloc+1:2*nloc, nloc+1:2*nloc) = zx42{k}(:,1:nloc,1:nloc);
  zx42_flipped{k}(:,nloc+1:2*nloc,1:nloc) = permute(zx42{k}(:,nloc+1:2*nloc,1:nloc), [1 3 2]);
  
  coh2_flipped{k}(:,1:nloc,1:nloc) = coh2{k}(:,nloc+1:2*nloc, nloc+1:2*nloc);
  coh2_flipped{k}(:,nloc+1:2*nloc, nloc+1:2*nloc) = coh2{k}(:,1:nloc,1:nloc);
  coh2_flipped{k}(:,nloc+1:2*nloc,1:nloc) = permute(coh2{k}(:,nloc+1:2*nloc,1:nloc), [1 3 2]);
  
  
  coh4_flipped{k}(:,1:nloc,1:nloc) = coh4{k}(:,nloc+1:2*nloc, nloc+1:2*nloc);
  coh4_flipped{k}(:,nloc+1:2*nloc, nloc+1:2*nloc) = coh4{k}(:,1:nloc,1:nloc);
  coh4_flipped{k}(:,nloc+1:2*nloc,1:nloc) = permute(coh4{k}(:,nloc+1:2*nloc,1:nloc), [1 3 2]);
  
  
  cohC{k} = (coh1{k} + coh4_flipped{k})./2;
  cohIC{k} = (coh3{k} + coh2_flipped{k})./2;
  
  zx{k} = (zx13{k} + zx42_flipped{k})./2;
end

allconnectivity = cell(n,1);
allcohC = cell(n,1);
allcohIC = cell(n,1);

for k=1:numel(zx)
  nloc = size(zx13{k},3)/2;
  x = ones(2*nloc);
  for subj=1:n
    tmpcoh = squeeze(zx{k}(subj,:,:));
    tmpcoh(find(triu(x))) = nan;
    tmpcoh = tmpcoh(:);
    tmpcoh = tmpcoh(~isnan(tmpcoh));
    allconnectivity{subj} = [allconnectivity{subj}, tmpcoh'];
    ncomparisson(k) = numel(tmpcoh);
    
    tmpC = squeeze(cohC{k}(subj,:,:));
    tmpC(find(triu(x))) = nan;
    tmpC = tmpC(:);
    tmpC = tmpC(~isnan(tmpC));
    allcohC{subj} = [allcohC{subj}, tmpC'];
    
    tmpIC = squeeze(cohIC{k}(subj,:,:));
    tmpIC(find(triu(x))) = nan;
    tmpIC = tmpIC(:);
    tmpIC = tmpIC(~isnan(tmpIC));
    allcohIC{subj} = [allcohIC{subj}, tmpIC'];
    
  end
end
allconnectivity = cat(1, allconnectivity{:});
allcohC = cat(1, allcohC{:});
allcohIC = cat(1, allcohIC{:});


data = [];
data.label = {'coh'};
data.time = 1:size(allconnectivity,2);
data.trial = reshape(allconnectivity, [n 1 size(allconnectivity,2)]);
data.dimord = 'rpt_chan_time';

nul = data;
nul.trial = nul.trial*0;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'trial';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';
stat = ft_timelockstatistics(cfgs, data, nul);

save(fullfile([alldir, sprintf('analysis/stat_coh_%s.mat', toi)]), 'stat', 'allconnectivity', 'data', 'ncomparisson', 'zx', 'within_vs_between_hemiC', 'zx_within', 'zx_between', 'within_CvsIC', 'between_CvsIC','coh_within_C', 'coh_within_IC', 'coh_between_C', 'coh_between_IC', 'allcohC', 'allcohIC');

