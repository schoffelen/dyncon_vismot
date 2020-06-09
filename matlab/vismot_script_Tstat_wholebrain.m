load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
if ~exist('stratifyflag', 'var'), stratifyflag = true; end
if ~exist('toi', 'var'), toi = 'post'; end

frequency = [6 10 22 38 42 58 62 78 82];
n=19;
dat1 = zeros(numel(sourcemodel.inside), numel(frequency), n);
dat2 = zeros(numel(sourcemodel.inside), numel(frequency), n);
raw = zeros(4,numel(sourcemodel.inside), numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    switch toi
      case 'post'
        if stratifyflag
          d = fullfile([datadir, sprintf('%s_source3d4mm_post_%03d_prewhitened_stratified.mat', list{m}, k)]);
        else
          d = fullfile([datadir, sprintf('%s_source3d4mm_post_%03d_prewhitened.mat', list{m}, k)]);
        end
      case 'pre'
        d = fullfile([datadir, sprintf('%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, k)]);
    end
    dum = load(d, 'stat', 'yuenD_13','yuenD_42');
    dum.stat = hemiflip(dum.stat, {'stat4', 'stat2'});
    dum.dim=dum.stat.dim;
    dum = hemiflip(dum, 'yuenD_42');
    raw(1,:,cnt+1,m) = dum.stat.stat1.*dum.yuenD_13;
    raw(2,:,cnt+1,m) = dum.stat.stat2.*dum.yuenD_42;
    raw(3,:,cnt+1,m) = dum.stat.stat3.*dum.yuenD_13;
    raw(4,:,cnt+1,m) = dum.stat.stat4.*dum.yuenD_42;
    dat1(:,cnt+1, m) = (dum.stat.stat1 + dum.stat.stat4)./2;
    dat2(:,cnt+1, m) = (dum.stat.stat3 + dum.stat.stat2)./2;
  end
  clear dum
  cnt=cnt+1;
end

% make data structure
freqs = { 'theta', 6, 6
  'alpha', 10, 10
  'beta', 22, 22
  'gamma1', [38 42], 40
  'gamma2', [58 62], 60
  'gamma3', [78 82], 80};

dat1 = permute(dat1, [3,1,2]);
dat2 = permute(dat2, [3,1,2]);
source1 = sourcemodel;
source1.dimord = 'rpt_pos_freq';
source1.freq = cat(2,freqs{:,end});
source2=source1;

% fill data structure. Average within bands.
m=0;
for k=1:size(freqs,1)
  freqidx = find(ismember(freqs{k,2}, frequency));
  lmin = m+1;
  m = max([m+freqidx]);
  freqidx = freqidx+(lmin-1);
  source1.stat(:,:,k) = nanmean(dat1(:,:,freqidx),3);
  source2.stat(:,:,k) = nanmean(dat2(:,:,freqidx),3);
end

% make the source's inside symmetrical
% make the source's inside symmetrical
tmp=zeros(source1.dim);
tmp(source1.inside)=1;
tmp=flip(tmp,1);
tmp(source1.inside)=tmp(source1.inside)+1;
tmp(tmp<2) = 0;
tmp(tmp==2) = 1;
source1.inside = tmp(:);
source1.stat(:,source1.inside==0,:) = nan;
source2.inside = tmp(:);
source2.stat(:,source2.inside==0,:) = nan;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'statfun_yuenTtest';
cfgs.yuen.type = 'depsamples';
cfgs.yuen.percent = 0.1;
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'cluster';
cfgs.numrandomization = 10000;
cfgs.clusteralpha = 0.05;
cfgs.correcttail = 'prob';
cfgs.clusterstatistic = 'maxsum';
cfgs.clusterthreshold = 'nonparametric_individual';
for k=1:6
  cfg=[];
  cfg.frequency = source1.freq(k);
  stat{k} = ft_sourcestatistics(cfgs, ft_selectdata(cfg, source1), ft_selectdata(cfg, source2));
end

% calculate raw effect size
tmpraw = nan(size(raw)); tmpraw(:,:,size(freqs,1)+1:end,:)=[];
tmpraw(:,:,1:3,:) = raw(:,:,1:3,:); % theta, alpha, beta
tmpraw(:,:,4,:) = nanmean(raw(:,:,4:5,:),3); % gamma 1
tmpraw(:,:,5,:) = nanmean(raw(:,:,6:7,:),3); % gamma 2
tmpraw(:,:,6,:) = nanmean(raw(:,:,8:9,:),3); % gamma 3
raw = tmpraw;

for k=1:numel(stat)
  maskpos{k} = find(stat{k}.mask==1 & stat{k}.stat>0);
  maskneg{k} = find(stat{k}.mask==1 & stat{k}.stat<0);
end

clear c ic
c = squeeze((raw(1,:,:,:) + raw(4,:,:,:))./2);
ic = squeeze((raw(3,:,:,:) + raw(2,:,:,:))./2);
c_ic = c./ic-1;

for k=1:6
  poseffectsize(:,k) = squeeze(mean(c_ic(maskpos{k},k,:),1));
  negeffectsize(:,k) = squeeze(mean(c_ic(maskneg{k},k,:),1));
end

for k=1:numel(stat)
  effectsize_largest_cluster_pos(:,k) = mean(c_ic(stat{k}.posclusterslabelmat==1,k,:),1);
  effectsize_largest_cluster_neg(:,k) = mean(c_ic(stat{k}.negclusterslabelmat==1,k,:),1);
end

switch toi
  case 'post'
    if stratifyflag
      filename = fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']);
    else
      filename = fullfile([alldir, 'analysis/stat_bf_post']);
    end
  case 'pre'
    filename = fullfile([alldir, 'analysis/stat_bf_pre_wholebrain']);
end
save(filename, 'cfgs', 'stat','c','ic', 'source1','source2', 'sourcemodel', 'freqs', 'poseffectsize', 'negeffectsize', 'effectsize_largest_cluster_pos','effectsize_largest_cluster_neg', 'c_ic')


