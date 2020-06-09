load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

sourceidx_orig = [55526, 57373];% 62382, 38502, 38637, 45730, 38524];
% low gam mot ipsi, contra, alpha par contra, occ ipsi, occ contra, mot
% ipsi, contra
            
x = find(sourcemodel.inside);
for k=1:numel(sourceidx_orig)
  sourceidx(k) = find(x==sourceidx_orig(k));
end

frequency = [38 42];
n=19;
cnt = 1;
cfg=[];
for k = 1:numel(frequency)
  for m = 1:n
        d = fullfile([datadir, sprintf('singletrialcorr/%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, frequency(k))]);
    load(d, 'freqpow');
    % hemiflip conditions 2 and 4
    freqpow.dim = sourcemodel.dim;
    idx = find(ismember(freqpow.trialinfo(:,4), [4 2]));
    tmpfreqpow = rmfield(freqpow, 'pow');
    tmpfreqpow.pow = zeros(numel(sourcemodel.inside), numel(idx));
    tmpfreqpow.pow(sourcemodel.inside,:) = freqpow.pow(:, idx);
    tmpfreqpow = hemiflip(tmpfreqpow, 'pow');
    % pool conditions again
    freqpow.pow(:, idx) = tmpfreqpow.pow(sourcemodel.inside,:);
    
    f{m}(:,:,k) = freqpow.pow(sourceidx,:);
    trlinfo{m} = freqpow.trialinfo;
  end
  cnt=cnt+1;
end

for k=1:n
  f{k} = mean(f{k},3);
end

clear idx
for k=1:n
    idx{1} = find((trlinfo{k}(:,1)==1 | trlinfo{k}(:,1)==3) & (trlinfo{k}(:,4)==1 |trlinfo{k}(:,4)==3));
    idx{2} = find((trlinfo{k}(:,1)==2 | trlinfo{k}(:,1)==4) & (trlinfo{k}(:,4)==2 |trlinfo{k}(:,4)==4));

  % correlate reaction times with power in the sensorimotor ROIs
  for ii=1:2
    r(k,ii,:) = corr(trlinfo{k}(idx{ii},3), (f{k}(:,idx{ii}))','type', 'spearman');
  end
end
% average over left and right hand responses
r = squeeze(mean(r,2));
[h p ci stats] = ttest(r)

save([alldir, 'analysis/stat_singletrialcorrelation'], 'stats', 'p', 'ci', 'h')
