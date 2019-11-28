load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
if ~exist('toi', 'var'), toi = 'post'; end
load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']))

cfg=[];
cfg.parameter = {'stat', 'mask'};
for k=1:6
    stat_int{k} = ft_sourceinterpolate(cfg, stat{k}, mri);
end

cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.maskparameter = 'mask2';
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
cfgp.opacitylim = [-0.1 0.1];
clim = [6 6 8 6 8 8];
for k=1:6
s=stat_int{k};
s.mask2 = s.stat;
xmin = min(s.stat); xmax = max(s.stat);
s.mask2(s.mask2<0.3*xmax & s.mask2>0.3*xmin & s.inside(:)==1)=0;
s.mask2(s.mask2>0)=1;
cfgp.funcolorlim = [-clim(k) clim(k)];
ft_sourceplot(cfgp, s);
end


% Note the FOIs and ROIs here:
%  ROIs in cm, for left hand response trials. For right hand response
%  trials, multiply first column with -1.
unit = 'cm';
roi = [    0.2   2.8   3.6   % theta frontal midline
          -2.2 -10.0  -0.4   % alpha occipital ipsi
           1.8  -9.2   0.4   % alpha occipital contra
          -3.8  -9.2   0.8   % beta occipital ipsi
          -1.4  -9.2   4.8   % beta parietal ipsi
          -4.2  -4.0   6.0   % beta premotor ipsi
           5.0  -2.4   5.6   % beta motor contra
          -3.4  -9.2   0.4   % gamma2 occipital ipsi
           3.0 -10.0   0.4   % gamma2 occipital contra
           3.0  -8.8   4.4   % gamma2 parietal contra
          -4.2  -5.2   6.8   % gamma2 motor ipsi
           6.2  -3.2   4.4   % gamma2 motor contra
          -2.2 -10.0   1.2   % gamma3 occipital ipsi
           3.0  -9.2   4.0   % gamma3 occipital contra
          -2.6  -7.6   6.4   % gamma3 parietal ipsi
           2.6  -6.0   8.0   % gamma3 parietal contra
           4.2  -3.2   6.8]; % gamma3 motor contra

description = [
{'theta frontal midline'}
{'alpha occipital ipsi'}
{'alpha occipital contra'}
{'beta occipital ipsi'}
{'beta parietal ipsi'}
{'beta premotor ipsi'}
{'beta motor contra'}
{'gamma2 occipital ipsi'}
{'gamma2 occipital contra'}
{'gamma2 parietal contra'}
{'gamma2 motor ipsi'}
{'gamma2 motor contra'}
{'gamma3 occipital ipsi'}
{'gamma3 occipital contra'}
{'gamma3 parietal ipsi'}
{'gamma3 parietal contra'}
{'gamma3 motor contra'}];
foi = [6 10 10 22 22 22 22 60 60 60 60 60 80 80 80 80 80]';

% create matrices that correspond to the connections for which we want to
% test coherence. That is: all connections where at least one of the two
% ROIs has an post-cue power effect (e.g. theta-ROI to all, but not
% all-to-all)
freqs = { 'theta', 6, 6
  'alpha', 10, 10
  'beta', 22, 22
  'gamma1', [38 42], 40
  'gamma2', [58 62], 60
  'gamma3', [78 82], 80};
  
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
    
    % do not test coherence between neighbouring freqs at nearby locations
    % (e.g. gamma1-gamma2 occ ipsi). Manually exclude:
    excludemanual = [2 4; 2 8; 2 13; 4 8; 4 13; 8 13; 3 9; 3 14; 9 14; 5 15; 6 11; 7 12; 7 17; 12 17; 10 16];
    for ii=1:size(excludemanual,1)
        roi_to_roi{k}(excludemanual(ii,1),excludemanual(ii,2)) = nan;
        roi_to_roi{k}(excludemanual(ii,2),excludemanual(ii,1)) = nan;
    end
end

save([alldir 'analysis/roi.mat'], 'roi', 'description', 'unit', 'foi', 'roi_to_roi');



% Find effect size in ROIs
% find ROI indices
idx = find_dipoleindex(sourcemodel, [roi; roi.*[-1 1 1]]);
freqs = [6 10 22 40 60 80];
foi = [foi; foi];
for k=1:numel(foi)
  freq_idx(k) = find(freqs==foi(k));
end
for m=1:numel(foi)
  effectsize_roi(:,m) = c_ic(idx(m), freq_idx(m),:);
end

for m=1:numel(freq_idx)
  tmpc(:,m) = c(idx(m), freq_idx(m), :);
  tmpic(:,m) = ic(idx(m), freq_idx(m), :);
end
ic_roi=tmpic;
c_roi=tmpc;

switch toi
  case 'post'
if stratifyflag
  filename = fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']);
else
  filename = fullfile([alldir, 'analysis/stat_bf_post.mat']);
end
  case 'pre'
    filename = fullfile([alldir, 'analysis/stat_bf_pre_wholebrain']);
end
save(filename, 'effectsize_roi', 'idx', 'freq_idx','c_roi', 'ic_roi', '-append')