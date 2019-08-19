load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
if ~exist('stratifyflag', 'var'), stratifyflag = false; end

% Pooled statistic: (left response C-IC minus right response C-IC)
whichstat = 'statResp';
frequency = [10 22 38 42 58 62 78 82];
n=19;
dat = zeros(numel(sourcemodel.inside), numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    if stratifyflag
      d = fullfile([datadir, sprintf('%s_source3d4mm_post_%03d_prewhitened_stratified.mat', list{m}, k)]);
    else
      d = fullfile([datadir, sprintf('%s_source3d4mm_post_%03d_prewhitened.mat', list{m}, k)]);
    end
    dum = load(d, 'stat');
    raw{m, cnt+1} = load(d,'source');
    dat(:,cnt+1, m) = dum.stat.(whichstat);
  end
  clear dum
  cnt=cnt+1;
end

% make data structure
foi = { 'alpha', 10, 10
        'beta', 22, 22
        'gamma1', [38 42], 40
        'gamma2', [58 62], 60
        'gamma3', [78 82], 80};

dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = cat(2,foi{:,end});

% fill data structure. Average within bands.
l=0;
for k=1:size(foi,1)
  freqidx = find(ismember(foi{k,2}, frequency));
  lmin = l+1;
  l = max([l+freqidx]);
  freqidx = freqidx+(lmin-1);
  source.stat(:,:,k) = nanmean(dat(:,:,freqidx),3);
end

% compare with zero
nul = source;
nul.stat(:)=0;
source.avg = squeeze(nanmean(source.stat,1));

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'cluster';
cfgs.numrandomization = 10000;
cfgs.clusteralpha = 0.025;
cfgs.correcttail = 'prob';
stat = ft_sourcestatistics(cfgs, source, nul);

% hemiflip raw 4-2 conditions
for k=1:n
  for f=1:8
    for c=[2 4]
      raw{k,f}.source(c).avg.dim = raw{k,f}.source(c).dim;
      raw{k,f}.source(c).avg = hemiflip(raw{k,f}.source(c).avg, 'pow');
    end
  end
end

% calculate raw effect size
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

mask = find(stat.mask==1);
posidx = find(stat.stat>0);
maskpos = mask(ismember(mask, posidx));
negidx = find(stat.stat<0);
maskneg = mask(ismember(mask, negidx));

clear c ic
for k=1:n
  c(k,:,:) = (raw{k,1}+raw{k,4})./2;
  ic(k,:,:) = (raw{k,2} + raw{k,3})./2;
end

for k=1:n
  c_ic(k,:,:) = (raw{k,1}./raw{k,3}-1 + raw{k,4}./raw{k,2}-1)./2;
end
c_ic2 = reshape(c_ic, n, []);

poseffectsize = mean(c_ic2(:,maskpos),2);
negeffectsize = mean(c_ic2(:,maskneg),2);

effectsize_largest_cluster = mean(c_ic(:,find(stat.posclusterslabelmat==1)),2);

% pool across hemispheres by subtracting the other hemisphere
stat_semhemi = stat;
for k=1:numel(stat.freq)
  tmpx=stat.stat(:,k);
  dum=nan(sourcemodel.dim);
  dum=reshape(tmpx,sourcemodel.dim);
  dum(~sourcemodel.inside)=nan;
  dum=dum-flip(dum,1);
  stat_semhemi.stat(:,k) = dum(:)/2;
  clear dum tmpx
end

if stratifyflag
  save(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']), 'stat', 'source', 'sourcemodel', 'foi', 'stat_semhemi', 'poseffectsize', 'negeffectsize', 'effectsize_largest_cluster', 'c_ic')
else
  save(fullfile([alldir, 'analysis/stat_bf_post.mat']), 'stat', 'source', 'sourcemodel', 'foi', 'stat_semhemi', 'poseffectsize', 'negeffectsize', 'effectsize_largest_cluster', 'c_ic')
end

%% Define ROI's by browsing through Ortho Maps
if ~stratifyflag
  cmap = flipud(brewermap(64,'RdBu'));
  cfgp=[];
  cfgp.funcolormap = cmap;
  cfgp.method     = 'ortho';
  cfgp.funparameter = 'stat';
  cfgp.location = 'max';
  for k=1:numel(stat_semhemi.freq)
    cfgp.frequency = stat.freq(k);
    ft_sourceplot(cfgp, stat_semhemi)
  end
end

% Note the FOIs and ROIs here:
roi = [
  {'REGION'}    {'FREQUENCY BAND'}  {'LOCATION LEFT'}   {'LOCATION RIGHT'}
  {'occipital'} {'alpha'}           {[-1.8 -9.6 0.4]}   {[1.8 -9.6 0.4]}
  {'occipital'} {'gamma1'}          {[-3.0 -10.4 0.0]}  {[3.0 -10.4 0.0]}
  {'occipital'} {'gamma2'}          {[-2.6 -9.2 0.8]}   {[2.6 -9.2 0.8]}
  {'occipital'} {'gamma3'}          {[-2.6 -10.0 0.8]}   {[2.6 -10.0 0.8]}
  {'parietal'}  {'alpha'}           {[-3.0 -7.6 5.6]}   {[3.0 -7.6 5.6]}
  {'parietal'}  {'gamma1'}          {[-3.4 -8.8 4.4]}   {[3.4 -8.8 4.4]}
  {'parietal'}  {'gamma2'}          {[-3.0 -7.6 6.0]}   {[3.0 -7.6 6.0]}
  {'parietal'}  {'gamma3'}          {[-2.6 -8.4 4.8]}   {[2.6 -8.4 4.8]}
  {'motor'}     {'alpha'}           {[-4.6 -0.4 6.4]}   {[4.6 -0.4 6.4]}
  {'motor'}     {'beta'}            {[-3.8 -4.0 5.6]}   {[3.8 -4.0 5.6]}
  {'motor'}     {'gamma1'}          {[-4.6 -3.6 5.6]}   {[4.6 -3.6 5.6]}
  {'motor'}     {'gamma3'}          {[-4.2 -3.6 7.2]}   {[4.2 -3.6 7.2]}];
foi = [
  {'FREQUENCY BAND'}  {'FREQUENCY RANGE'} {'CENTER FREQUENCY'} {'REPRESENTATIVE FREQUENCIES'} {'SMOOTHING'}
  {'alpha' }          {[8 12]}            {10}                 {10}                           {2}
  {'beta'  }          {[14 30]}           {22}                 {22}                           {8}
  {'gamma1'}          {[30 50]}           {40}                 {[38 42]}                      {8}
  {'gamma2'}          {[50 70]}           {60}                 {[58 62]}                      {8}
  {'gamma3'}          {[70 90]}           {80}                 {[78 82]}                      {8}];

unit = 'cm';
if ~stratifyflag
  save('/project/3011085.03/analysis/source/roi.mat', 'roi', 'foi', 'unit');
end

% Find effect size in ROIs
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


for m=1:numel(freq_idx)
  effectsize_roi_left(:,m) = c_ic(:,l(m), freq_idx(m));
  effectsize_roi_right(:,m) =  c_ic(:,r(m), freq_idx(m));
end

for m=1:numel(freq_idx)
  c_left(:,m) = c(:,l(m), freq_idx(m));
  ic_left(:,m) = ic(:,l(m), freq_idx(m));
  c_right(:,m) = c(:,r(m), freq_idx(m));
  ic_right(:,m) = ic(:,r(m), freq_idx(m));
end

if stratifyflag
  save(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']), 'effectsize_roi_right','effectsize_roi_left', 'l', 'r', 'freq_idx','c', 'ic', 'c_left', 'c_right', 'ic_left', 'ic_right', '-append')
else
  save(fullfile([alldir, 'analysis/stat_bf_post.mat']), 'effectsize_roi_right','effectsize_roi_left', 'l', 'r', 'freq_idx', 'freq_idx','c', 'ic', 'c_left', 'c_right', 'ic_left', 'ic_right', '-append')
end
