load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

% Pooled statistic: (left response C-IC minus right response C-IC)
whichstat = 'statResp';
frequency = [10 22 38 42 58 62 78 82];%[10 16 24 26 48:4:92];%[4:2:30 40:4:100];
n=19;
dat = zeros(numel(sourcemodel.inside), numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%s_source3d4mm_post_%03d_prewhitened.mat', list{m}, k)]);
    dum = load(d, 'stat');
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
cfgs.numrandomization = 1000;
cfgs.clusteralpha = 0.025;
cfgs.correcttail = 'prob';
stat = ft_sourcestatistics(cfgs, source, nul);

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
save(fullfile([alldir, 'analysis/stat_bf_post.mat']), 'stat', 'source', 'sourcemodel', 'foi', 'stat_semhemi')


%% Define ROI's by browsing through Ortho Maps
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

% Note the FOIs and ROIs here:
roi = [
  {'REGION'}    {'FREQUENCY BAND'}  {'LOCATION LEFT'}   {'LOCATION RIGHT'}
  {'occipital'} {'alpha'}           {[-2.6 -9.6 3.2]}   {[2.6 -9.6 3.2]}
  {'occipital'} {'gamma1'}          {[-3.0 -10.4 0.8]}  {[3.0 -10.4 0.8]}
  {'occipital'} {'gamma2'}          {[-2.6 -9.2 0.8]}   {[2.6 -9.2 0.8]}
  {'occipital'} {'gamma3'}          {[-3.4 -9.2 1.6]}   {[3.4 -9.2 1.6]}
  {'parietal'}  {'alpha'}           {[-3.4 -7.6 5.6]}   {[3.4 -7.6 5.6]}
  {'parietal'}  {'gamma1'}          {[-3.4 -8.8 4.4]}   {[3.4 -8.8 4.4]}
  {'parietal'}  {'gamma2'}          {[-3.4 -7.6 5.6]}   {[3.4 -7.6 5.6]}
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

save('/project/3011085.03/analysis/source/roi.mat', 'roi', 'foi');


