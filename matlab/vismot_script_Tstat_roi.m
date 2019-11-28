load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

frequency = [6 10 22 38 42 58 62 78 82];
n=19;
dat1 = zeros(numel(sourcemodel.inside), numel(frequency), n);
dat2 = zeros(numel(sourcemodel.inside), numel(frequency), n);
raw = zeros(4,numel(sourcemodel.inside), numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, k)]);
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

m=0;
for k=1:size(freqs,1)
  freqidx = find(ismember(freqs{k,2}, frequency));
  lmin = m+1;
  m = max([m+freqidx]);
  freqidx = freqidx+(lmin-1);
  source1.stat(:,:,k) = nanmean(dat1(:,:,freqidx),3);
  source2.stat(:,:,k) = nanmean(dat2(:,:,freqidx),3);
end


clear foi
load(fullfile([alldir, 'analysis/roi.mat']));

roi_idx = find_dipoleindex(sourcemodel, roi);
for k=1:size(roi,1)
  tmpidx = find(foi(k)==cat(1,freqs{:, end}));
  tmpstat_roiC(:,k) = source1.stat(:,roi_idx(k),tmpidx);
  tmpstat_roiIC(:,k) = source2.stat(:,roi_idx(k),tmpidx);
end

s1 =[];
s1.stat(:,1,:) = tmpstat_roiC;
s1.time = 1:size(tmpstat_roiC,2);
s1.dimord = 'rpt_chan_time';
s1.label{1} = 'pre_roi_pow';

s2 = rmfield(s1, 'stat');
s2.stat(:,1,:) = tmpstat_roiIC;

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
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';

stat = ft_timelockstatistics(cfgs, s1, s2);

% calculate raw effect size
tmpraw = nan(size(raw)); tmpraw(:,:,size(freqs,1)+1:end,:)=[];
tmpraw(:,:,1:3,:) = raw(:,:,1:3,:); % theta, alpha, beta
tmpraw(:,:,4,:) = nanmean(raw(:,:,4:5,:),3); % gamma 1
tmpraw(:,:,5,:) = nanmean(raw(:,:,6:7,:),3); % gamma 2
tmpraw(:,:,6,:) = nanmean(raw(:,:,8:9,:),3); % gamma 3
raw = tmpraw;
clear tmpraw

tmpraw = nan(4, n, size(roi,1));
for k=1:size(roi,1)
  tmpidx = find(foi(k)==cat(1,freqs{:, end}));
  tmpraw(:,:,k) = squeeze(raw(:,roi_idx(k),tmpidx,:));
end
raw=tmpraw;

c = squeeze((raw(1,:,:) + raw(4,:,:))./2);
ic = squeeze((raw(3,:,:) + raw(2,:,:))./2);
c_ic = c./ic-1;

filename = [alldir, 'analysis/stat_bf_pre.mat'];
save(filename, 'c', 'ic', 'c_ic','description', 'sourcemodel','stat', 'stat_roi', 'd', 'raw');


