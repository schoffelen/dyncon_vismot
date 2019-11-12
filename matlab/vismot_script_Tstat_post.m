load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
if ~exist('stratifyflag', 'var'), stratifyflag = true; end

% Pooled statistic: (left response C-IC minus right response C-IC)
whichstat = 'statResp';
frequency = [6 10 22 38 42 58 62 78 82];
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
foi = { 'theta', 6, 6
  'alpha', 10, 10
  'beta', 22, 22
  'gamma1', [38 42], 40
  'gamma2', [58 62], 60
  'gamma3', [78 82], 80};

dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = cat(2,foi{:,end});

% fill data structure. Average within bands.
m=0;
for k=1:size(foi,1)
  freqidx = find(ismember(foi{k,2}, frequency));
  lmin = m+1;
  m = max([m+freqidx]);
  freqidx = freqidx+(lmin-1);
  source.stat(:,:,k) = nanmean(dat(:,:,freqidx),3);
end

% compare with zero
nul = source;
nul.stat = nul.stat*0;
source.avg = squeeze(nanmean(source.stat,1));

% make the source's inside symmetrical
tmp=zeros(source.dim);
tmp(source.inside)=1;
tmp=flip(tmp,1);
tmp(source.inside)=tmp(source.inside)+1;
for k=1:19
  for m = 1:6
    tmp=reshape(source.stat(k,:,m),source.dim);
    source.stat(k,:,m)=tmp(:);
  end
end


cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'statfun_yuenTtest';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'cluster';
cfgs.numrandomization = 10000;
cfgs.clusteralpha = 0.025;
cfgs.correcttail = 'prob';
cfgs.clusterstatistic = 'max';
cfgs.clusterthreshold = 'nonparametric_individual';
for k=1:6
  cfg=[];
  cfg.frequency = source.freq(k);
  stat{k} = ft_sourcestatistics(cfgs, ft_selectdata(cfg, source), ft_selectdata(cfg, nul));
end

% hemiflip raw 4-2 conditions
for k=1:n
  for f=1:9
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
    tmpraw{k,c}(:,2) = raw{k,2}.source(c).avg.pow; % alpha
    tmpraw{k,c}(:,3) = raw{k,3}.source(c).avg.pow; % beta
    tmpraw{k,c}(:,4) = nanmean([raw{k,4}.source(c).avg.pow, raw{k,5}.source(c).avg.pow],2); % gamma 1
    tmpraw{k,c}(:,5) = nanmean([raw{k,6}.source(c).avg.pow, raw{k,7}.source(c).avg.pow],2); % gamma 2
    tmpraw{k,c}(:,6) = nanmean([raw{k,8}.source(c).avg.pow, raw{k,9}.source(c).avg.pow],2); % gamma 3
  end
end
raw = tmpraw;

for k=1:numel(stat)
  maskpos{k} = find(stat{k}.mask==1 & stat{k}.stat>0);
  maskneg{k} = find(stat{k}.mask==1 & stat{k}.stat<0);
end

clear c ic
for k=1:n
  c(k,:,:) = (raw{k,1}+raw{k,4})./2;
  ic(k,:,:) = (raw{k,2} + raw{k,3})./2;
end

for k=1:n
  c_ic(k,:,:) = (raw{k,1}./raw{k,3}-1 + raw{k,4}./raw{k,2}-1)./2;
end

for k=1:6
  poseffectsize(:,k) = mean(c_ic(:,maskpos{k},k),2);
  negeffectsize(:,k) = mean(c_ic(:,maskneg{k},k),2);
end

for k=1:numel(stat)
  effectsize_largest_cluster_pos(:,k) = mean(c_ic(:,(stat{k}.posclusterslabelmat==1),k),2);
  effectsize_largest_cluster_neg(:,k) = mean(c_ic(:,(stat{k}.negclusterslabelmat==1),k),2);
end

% pool across hemispheres by subtracting the other hemisphere
for k=1:numel(stat)
  stat_semhemi{k} = stat{k};
  
  tmpx=stat{k}.stat;
  dum=nan(sourcemodel.dim);
  dum=reshape(tmpx,sourcemodel.dim);
  dum(~sourcemodel.inside)=nan;
  dum=dum-flip(dum,1);
  stat_semhemi{k}.stat = dum(:)/2;
  clear dum tmpx
end

if stratifyflag
  save(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']), 'stat','c','ic', 'source', 'sourcemodel', 'foi', 'stat_semhemi', 'poseffectsize', 'negeffectsize', 'effectsize_largest_cluster_pos','effectsize_largest_cluster_neg', 'c_ic')
else
  save(fullfile([alldir, 'analysis/stat_bf_post.mat']), 'stat', 'source', 'sourcemodel', 'foi', 'stat_semhemi', 'poseffectsize', 'negeffectsize', 'effectsize_largest_cluster_pos','effectsize_largest_cluster_neg', 'c_ic')
end

%% Define ROI's by browsing through Ortho Maps
if stratifyflag
  cmap = flipud(brewermap(64,'RdBu'));
  cfgp=[];
  cfgp.funcolormap = cmap;
  cfgp.method     = 'ortho';
  cfgp.funparameter = 'stat';
  cfgp.location = 'max';
  for k=1:numel(stat)
    ft_sourceplot(cfgp, stat_semhemi{k})
  end
end

% Note the FOIs and ROIs here:
% roi for left hand response, unit cm
roi = [-0.6   2.0  3.2 % theta frontal midline
  1.0  -8.8  0.4 % alpha occipital contra
  -1.4  -9.6  2.4 % alpha occipital ipsi
  1.4  -4.4  8.4 % alpha parietal contra
  -4.2  -1.6  6.4 % alpha motor ipsi
  -3.4  -8.8  4.0 % beta occipital ipsi
  -5.8   1.6  4.0 % beta premotor ipsi
  3.4  -2.8  8.0 % beta motor contra
  -4.2  -4.4  5.6 % beta motor ipsi
  -2.6  -8.8  0.4 % gamma1 occipital ipsi
  4.2  -7.6  4.0 % gamma1 par contra
  2.6  -10.0 0.0 % gamma2 occipital contra
  -3.8  -9.6 -0.4 % gamma2 occipital ipsi
  3.4  -7.2  5.6 % gamma2 parietal contra
  -3.8  -5.6  6.8 % gamma2 parietal ipsi
  3.0  -5.6  8.0 % gamma3 parietal contra
  4.6  -2.4  6.4 % gamma3 motor contra
  -2.6  -2.4  7.6 ];% gamma3 motor ipsi

description = [{'theta frontal midline'}
  {'alpha occipital contra'}
  {'alpha occipital ipsi'}
  {'alpha parietal contra'}
  {'alpha motor ipsi'}
  {'beta occipital ipsi'}
  {'beta premotor ipsi'}
  {'beta motor contra'}
  {'beta motor ipsi'}
  {'gamma1 occipital ipsi'}
  {'gamma1 par contra'}
  {'gamma2 occipital contra'}
  {'gamma2 occipital ipsi'}
  {'gamma2 parietal contra'}
  {'gamma2 parietal ipsi'}
  {'gamma3 parietal contra'}
  {'gamma3 motor contra'}
  {'gamma3 motor ipsi'}];
foi = [6 10 10 10 10 22 22 22 22 40 40 60 60 60 60 80 80 80]';

% create matrices that correspond to the connections for which we want to
% test coherence. That is: all connections where at least one of the two
% ROIs has an post-cue power effect (e.g. theta-ROI to all, but not
% all-to-all)
for k=1:size(freqs,1)
  roi_to_roi{k} = nan(size(roi,1));
end
f = cat(1,freqs{:,3});
for k=1:size(freqs,1)
  tmpidx = find(foi==f(k));
  for m=1:numel(tmpidx)
    roi_to_roi{k}(tmpidx(m),:)=1;
    roi_to_roi{k}(:,tmpidx(m))=1;
  end
  % do not take the diagonal
  dum = diag(diag(nan(size(roi,1)))); dum(~isnan(dum))=1;
  roi_to_roi{k} = roi_to_roi{k}.*dum;
  
  % take only the lower part of the triangle
  dum = ones(size(roi,1));
  roi_to_roi{k}(find(triu(dum))) = nan;
end
unit = 'cm';
if stratifyflag
  save('/project/3011085.03/analysis/roi.mat', 'roi', 'description', 'unit', 'foi', 'roi_to_roi');
end


% Find effect size in ROIs
% find ROI indices
idx = find_dipoleindex(sourcemodel, [roi; roi.*[-1 1 1]]);
freqs = [6 10 22 40 60 80];
foi = [foi; foi];
for k=1:numel(foi)
  freq_idx(k) = find(freqs==foi(k));
end
for m=1:numel(foi)
  effectsize_roi(:,m) = c_ic(:,idx(m), freq_idx(m));
end

for m=1:numel(freq_idx)
  tmpc(:,m) = c(:,idx(m), freq_idx(m));
  tmpic(:,m) = ic(:,idx(m), freq_idx(m));
end
ic=tmpic;
c=tmpc;

if stratifyflag
  filename = fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']);
else
  filename = fullfile([alldir, 'analysis/stat_bf_post.mat']);
end
save(filename, 'effectsize_roi', 'idx', 'freq_idx','c', 'ic', '-append')