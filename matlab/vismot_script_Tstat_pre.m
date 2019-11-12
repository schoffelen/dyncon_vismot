load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = {'stat13', 'stat42'};
frequency = [6 10 22 38 42 58 62 78 82];
n=19;
dat{1} = zeros(74784, numel(frequency), n);
dat{2} = zeros(74784, numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, k)]);
    dum = load(d, 'stat');
    raw{m, cnt+1} = load(d,'source');
    for l=1:2
      dat{l}(:, cnt+1,m) = dum.stat.(whichstat{l});
    end
  end
  clear dum
  cnt=cnt+1;
end
freqs = { 'theta', 6, 6
  'alpha', 10, 10
  'beta', 22, 22
  'gamma1', [38 42], 40
  'gamma2', [58 62], 60
  'gamma3', [78 82], 80};

for l=1:numel(dat)
  dat{l} = permute(dat{l}, [3,1,2]);
  source{l} = sourcemodel;
  source{l}.dimord = 'rpt_pos_freq';
  source{l}.freq = cat(2,freqs{:, end});
  
  m=0;
  for k=1:size(freqs,1)
    freqidx = find(ismember(freqs{k,2}, frequency));
    lmin = m+1;
    m = max([m+freqidx]);
    freqidx = freqidx+(lmin-1);
    source{l}.stat(:,:,k) = nanmean(dat{l}(:,:,freqidx),3);
  end
end
clear dat

clear foi
load(fullfile([alldir, 'analysis/roi.mat']));

roi_idx = find_dipoleindex(sourcemodel, [roi; roi.*[-1 1 1]]);
for k=1:size(roi,1)
  tmpidx = find(foi(k)==cat(1,freqs{:, end}));
  tmpstat_roi13(:,k) = source{1}.stat(:,roi_idx(k),tmpidx);
  tmpstat_roi42(:,k) = source{2}.stat(:,roi_idx(k+size(roi,1)),tmpidx);
end
stat_roi = (tmpstat_roi13 + tmpstat_roi42)./2;

d =[];
d.stat(:,1,:) = stat_roi;
d.time = 1:size(stat_roi,2);
d.dimord = 'rpt_chan_time';
d.label{1} = 'pre_roi_pow';

nul = d;
nul.stat=0*nul.stat;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'statfun_yuenTtest';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';

stat = ft_timelockstatistics(cfgs, d, nul);

% average raw power within frequency bands for each condition
for k=1:n
  for c=1:4
    tmpraw{k,c}(:,1) = raw{k,1}.source(c).avg.pow; % theta
    tmpraw{k,c}(:,2) = raw{k,2}.source(c).avg.pow; % alpha
    tmpraw{k,c}(:,3) = raw{k,3}.source(c).avg.pow; % beta
    tmpraw{k,c}(:,4) = nanmean([raw{k,4}.source(c).avg.pow, raw{k,5}.source(c).avg.pow],2); % gamma 1
    tmpraw{k,c}(:,5) = nanmean([raw{k,6}.source(c).avg.pow, raw{k,7}.source(c).avg.pow],2); % gamma 2
    tmpraw{k,c}(:,6) = nanmean([raw{k,8}.source(c).avg.pow, raw{k,9}.source(c).avg.pow],2); % gamma 3
  end
end
raw = tmpraw;
clear tmpraw

tmpraw = nan(4, n, size(roi,1));
for s=1:n
  for k=1:size(roi,1)
    tmpidx = find(foi(k)==cat(1,freqs{:, end}));
    for l=[1 3]
      tmpraw(l,s,k) = raw{s,l}(roi_idx(k),tmpidx);
    end
    for l=[2 4]
      tmpraw(l,s,k) = raw{s,l}(roi_idx(k+size(roi,1)),tmpidx);
    end
  end
end


clear c ic

c = (squeeze(tmpraw(1,:,:)+tmpraw(4,:,:)))./2;
ic = (squeeze(tmpraw(3,:,:)+tmpraw(2,:,:)))./2;
c_ic = c./ic-1;

save('/project/3011085.03/analysis/stat_bf_pre.mat', 'c', 'ic', 'c_ic','description', 'sourcemodel','stat', 'stat_roi', 'd');


