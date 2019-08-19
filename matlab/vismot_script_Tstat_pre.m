load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = 'statResp';
frequency = [10 22 38 42 58 62 78 82];
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 0;
for k = fliplr(frequency)
  for m = 1:n
    d = fullfile([datadir, sprintf('%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, k)]);
    dum = load(d, 'stat');
    raw{m, cnt+1} = load(d,'source');
    dat(:, cnt+1,m) = dum.stat.(whichstat);
  end
  clear dum
  cnt=cnt+1;
end
foi = { 'alpha', 10, 10
        'beta', 22, 22
        'gamma1', [38 42], 40
        'gamma2', [58 62], 60
        'gamma3', [78 82], 80};

dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = cat(2,foi{:,end});
l=0;
for k=1:size(foi,1)
  freqidx = find(ismember(foi{k,2}, frequency));
  lmin = l+1;
  l = max([l+freqidx]);
  freqidx = freqidx+(lmin-1);
  source.stat(:,:,k) = nanmean(dat(:,:,freqidx),3);
end


load(fullfile([alldir, 'analysis/source/roi.mat']));
l = zeros(size(roi,1)-1,3);
r = zeros(size(roi,1)-1,3);
% find ROI indices
for k=1:size(roi,1)-1
  l(k,:) = roi{k+1,3};
  r(k,:) = roi{k+1,4};
  freq_idx(k) = find(strcmp(roi{k+1,2}, foi([2:end],1)));
end
l = find_dipoleindex(sourcemodel, l);
r = find_dipoleindex(sourcemodel, r);

stat_roi = nan(n,numel(l), 1);
for k=1:numel(l)
  lpow(:,k) = source.stat(:,l(k),freq_idx(k));
  rpow(:,k) = source.stat(:,r(k),freq_idx(k));
end
stat_roi = [lpow rpow];

d =[];
d.stat(:,1,:) = stat_roi;
d.time = 1:size(stat_roi,2);
d.dimord = 'rpt_chan_time';
d.label{1} = 'pre_roi_pow';

nul = d;
nul.stat=0*nul.stat;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';
stat = ft_timelockstatistics(cfgs, d, nul);

% hemiflip raw 4-2 conditions
for k=1:n
  for f=1:8
    for c=[2 4]
      raw{k,f}.source(c).avg.dim = raw{k,f}.source(c).dim;
      raw{k,f}.source(c).avg = hemiflip(raw{k,f}.source(c).avg, 'pow');
    end
  end
end
% average within frequency bands
for k=1:n
  for c=1:4
    tmpraw{k,c}(:,1) = raw{k,1}.source(c).avg.pow; % alpha
    tmpraw{k,c}(:,2) = raw{k,2}.source(c).avg.pow; % beta
    tmpraw{k,c}(:,3) = nanmean([raw{k,3}.source(c).avg.pow, raw{k,4}.source(c).avg.pow],2); % gamma 1
    tmpraw{k,c}(:,4) = nanmean([raw{k,5}.source(c).avg.pow, raw{k,6}.source(c).avg.pow],2); % gamma 2
    tmpraw{k,c}(:,5) = nanmean([raw{k,7}.source(c).avg.pow, raw{k,8}.source(c).avg.pow],2); % gamma 3
  end  
end
raw = tmpraw;
clear c ic
for k=1:n
  c(k,:,:) = (raw{k,1}+raw{k,4})./2;
  ic(k,:,:) = (raw{k,2} + raw{k,3})./2;
end

for m=1:numel(freq_idx)
  c_left(:,m) = c(:,l(m), freq_idx(m));
  ic_left(:,m) = ic(:,l(m), freq_idx(m));
  c_right(:,m) = c(:,r(m), freq_idx(m));
  ic_right(:,m) = ic(:,r(m), freq_idx(m));
end

for k=1:n
  c_ic(k,:,:) = (raw{k,1}+raw{k,4})./(raw{k,3}+raw{k,2})-1;
end
for m=1:numel(freq_idx)
  effectsize_roi_left(:,m) = c_ic(:,l(m), freq_idx(m));
  effectsize_roi_right(:,m) =  c_ic(:,r(m), freq_idx(m));
end

save('/project/3011085.03/analysis/stat_bf_pre.mat','c_left', 'ic_left', 'c_right', 'ic_right', 'c', 'ic', 'c_ic', 'l','r','freq_idx', 'sourcemodel', 'source','stat', 'stat_roi', 'd', 'rpow', 'lpow','effectsize_roi_left', 'effectsize_roi_right');


