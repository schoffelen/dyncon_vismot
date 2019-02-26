load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

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

save(fullfile([alldir, 'analysis/stat_bf_post.mat']), 'stat', 'source', 'sourcemodel', 'foi')

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


%% Define ROI's by browsing through Ortho Maps
source = ft_convert_units(source, 'mm');

cfgp=[];
cfgp.funcolormap = cmap;
cfgp.maskstyle = 'colormix';
cfgp.method     = 'ortho';
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'stat';

cfgp.frequency = frequency(1);
cfgp.location='max';
ft_sourceplot(cfgp, source);


%% Make tval bar graph of defined ROIs and FOIs
load('/project/3011085.03/analysis/source/roi.mat', 'ROI');
filename = fullfile(['/project/3011085.03/', 'analysis/', 'freq/', 'post_cue_pow_stat.mat']);
load(filename, 'source_int', 'source', 'stat', 'frequency','foi');

for k=1:3
    idx_left(k) = find_dipoleindex(source, ROI{k,2});
    idx_right(k) = find_dipoleindex(source, ROI{k,3});
end


figure;
% occipital
subplot(1,2,1);
y = [source.stat(idx_left(1),2), source.stat(idx_left(1),5); source.stat(idx_right(1),2), source.stat(idx_right(1),5)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1'}, 'location', 'northwest')
title('40-60 Hz')

% parietal
subplot(1,3,2);
y = [source.stat(idx_left(2),2), source.stat(idx_left(2),5) source.stat(idx_left(2),6); source.stat(idx_right(2),2), source.stat(idx_right(2),5) source.stat(idx_right(2),6)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1', 'gamma2'}, 'location', 'northwest')
title('parietal')

% motor
subplot(1,3,3);
y = [source.stat(idx_left(3),4), source.stat(idx_left(3),6); source.stat(idx_right(3),4), source.stat(idx_right(3),6)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'beta2', 'gamma2'}, 'location', 'northwest')
title('motor')
