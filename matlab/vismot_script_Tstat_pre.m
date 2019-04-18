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
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%s_source3d4mm_pre_%03d_prewhitened.mat', list{m}, k)]);
    dum = load(d, 'stat');
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

l = [
-2.6 -9.6 3.2	% occ alpha
-3.0 -10.4 0.8 % occ gam1
-2.6 -9.2 0.8 % occ gam2
-3.4 -9.2 1.6 % occ gam3
-3.4 -7.6 5.6 % par alpha
-3.4 -8.8 4.4 % par gam1
-3.4 -7.6 5.6 % par gam2
-2.6 -8.4 4.8 % par gam3
-4.6 -0.4 6.4 % mot alpha
-3.8 -4.0 5.6 % mot beta 
-4.6 -3.6 5.6 % mot gam1
-4.2 -3.6 7.2]; % mot gam3

freq_idx = [1 3 4 5 1 3 4 5 1 2 3 5]';

r=l; 
r(:,1) = -r(:,1);
l = find_dipoleindex(sourcemodel, l);
r = find_dipoleindex(sourcemodel, r);

stat_roi = nan(n,numel(l), 1);
for k=1:numel(l)
  stat_roi(:,k) = (source.stat(:,l(k),freq_idx(k))-source.stat(:,r(k),freq_idx(k)))./2;
end

source.freq=0;
source.stat(:,:,2:end)=[];

source.stat(:)=nan;
source.stat(:,1:numel(l),:) = stat_roi;
source.inside(:) = 0;
source.inside(1:numel(l))=1;

nul = source;
nul.stat=0*nul.stat;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferoni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';
stat = ft_sourcestatistics(cfgs, source, nul);

save('/project/3011085.03/analysis/stat_bf_pre.mat', 'l','r','freq_idx', 'sourcemodel', 'source','stat');


