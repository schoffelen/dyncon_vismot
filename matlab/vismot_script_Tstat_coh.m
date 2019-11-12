
load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
frequency = [6 10 22 38 42 58 62 78 82];
n=numel(list);

toi = 'pre';
roi_to = 'roi';
load(fullfile([alldir, 'analysis/roi.mat']));
nregions = 2*(size(roi,1));

freqs = { 'theta', 6, 6
  'alpha', 10, 10
  'beta', 22, 22
  'gamma1', [38 42], 40
  'gamma2', [58 62], 60
  'gamma3', [78 82], 80};

% load in data and average within frequency bands
cnt=1;
for k=1:6
  tmpfreqs = freqs{k,2};
  sz = numel(tmpfreqs);
  for subj=1:n
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

 % ROI  indices are dependent on response hand. Take the first half of 
 % indices for conditions 1/3, the second half for conditons 4/2. 
numroi = size(roi,1);
for k=1:size(freqs,1)
  zx13{k} = zx13{k}(:,1:numroi,1:numroi);
  zx42{k} = zx42{k}(:,numroi+1:2*numroi, numroi+1:2*numroi);
  coh1{k} = coh1{k}(:,1:numroi,1:numroi);
  coh2{k} = coh2{k}(:,numroi+1:2*numroi, numroi+1:2*numroi);
  coh3{k} = coh3{k}(:,1:numroi,1:numroi);
  coh4{k} = coh4{k}(:,numroi+1:2*numroi, numroi+1:2*numroi);
end

%% Is there any connectivity pattern significant individually (i.e. larger/smaller for C than IC)?
% combine 1-3 and 4-2
for k=1:size(freqs,1)  
  cohC{k} = (coh1{k} + coh4{k})./2; % raw coherence
  cohIC{k} = (coh3{k} + coh2{k})./2;
  
  zx{k} = (zx13{k} + zx42{k})./2; % stat
end

allconnectivity = cell(n,1);
allcohC = cell(n,1);
allcohIC = cell(n,1);

% take only the connections we're interested in (at least one ROI should
% have an effect in post cue window in that frequency).
for k=1:numel(zx)
  sel = roi_to_roi{k};
  
  for subj=1:n
    tmpcoh = squeeze(zx{k}(subj,:,:));
    tmpcoh = tmpcoh(sel>0);
    allconnectivity{subj} = [allconnectivity{subj}, tmpcoh'];
    ncomparisson(k) = sum(sel(:)>0);
    
    tmpC = squeeze(cohC{k}(subj,:,:));
    tmpC = tmpC(sel>0);
    allcohC{subj} = [allcohC{subj}, tmpC'];
    
    tmpIC = squeeze(cohIC{k}(subj,:,:));
    tmpIC = tmpIC(sel>0);
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
cfgs.statistic = 'statfun_yuenTtest';
cfgs.parameter = 'trial';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 20000;
cfgs.correcttail = 'prob';
stat = ft_timelockstatistics(cfgs, data, nul);

% structure the results into summary statistics.
if 1%any(stat.prob<cfgs.alpha)
  x = 1:numel(stat.stat);%find(stat.prob<cfgs.alpha);
  cmsm = cumsum(ncomparisson);
  dum = reshape(1:size(roi,1).^2, size(roi,1), size(roi,1));
  for k=1:numel(x)
    tmpfreqidx = find(x(k)<=cmsm, 1, 'first');
    effect(k).freq = freqs{tmpfreqidx,1};
    if tmpfreqidx==1
      selidx = x(k);
    else
      selidx = x(k)-cmsm(tmpfreqidx-1); 
    end
    idx_allconn = find(roi_to_roi{tmpfreqidx}>0);
    selidx = idx_allconn(selidx); % this is the index in the roi_to_roi that indicates the connection
    
    % now find the row and column
    [row, col] = find(dum==selidx);
    effect(k).roi1 = description{row};
    effect(k).roi2 = description{col};
    effect(k).stat = stat.stat(x(k));
    effect(k).prob_uncorrected = stat.prob(x(k));
    effect(k).ncomparisson = sum(ncomparisson);
    effect(k).correction = 'bonferroni';
    effect(k).prob_corrected = sum(ncomparisson)*stat.prob(x(k));
    effect(k).rawcoh_C = allcohC(:,x(k));
    effect(k).rawcoh_IC = allcohIC(:,x(k)); 
  end
  idx_sign_uncorrected = find(cat(1,effect(:).prob_uncorrected)<0.05);
  idx_sign_corrected = find(cat(1,effect(:).prob_corrected)<0.05);
else
  effect=[];
end
    
save(fullfile([alldir, sprintf('analysis/stat_coh_%s.mat', toi)]), 'stat', 'allconnectivity', 'data', 'ncomparisson', 'zx', 'allcohC', 'allcohIC', 'effect', 'idx_sign_corrected', 'idx_sign_uncorrected');

