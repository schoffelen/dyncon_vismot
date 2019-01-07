load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = 'statResp';
frequency = [6 10 16 24 26 48 52 58 62 78 82];%[10 16 24 26 48:4:92];%[4:2:30 40:4:100];
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%ssource3d4mm_post_%03d.mat', list{m}, k)]);
    dum = load(d, whichstat);
    tmp(m) = dum.(whichstat);
  end
  dat(:, cnt+1,1:n) = cat(2,tmp.stat);
  clear tmp dum
  cnt=cnt+1;
end
foi = {'theta', 6, 6
        'alpha', 10, 10
        'beta1', 16, 16
        'beta2', [24 26], 25
        'gamma1', [58 62], 50
        'gamma2', [78 82], 60
        'gamma3', [48 52], 80};
frequency = [6 10 16 25 50 60 80];
    
dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = frequency;
source.stat(:,:,1:3) = dat(:,:,1:3);
source.stat(:,:,4) = nanmean(dat(:,:,4:5),3);
source.stat(:,:,5) = nanmean(dat(:,:,6:7),3);
source.stat(:,:,6) = nanmean(dat(:,:,8:9),3);
source.stat(:,:,7) = nanmean(dat(:,:,10:11),3);
nul = source;
nul.stat(:)=0;
source_avg = rmfield(source, {'stat'});
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
cfgs.clusteralpha = 0.01;
cfgs.correcttail = 'no';
for k=1:size(foi,1)
    cfgs.frequency = source.freq(k);
    stat{k} = ft_sourcestatistics(cfgs, source, nul);
end

source = rmfield(source, 'stat');
for k=1:7
    source.stat(:,k) = stat{k}.stat;
    source.mask(:,k) = stat{k}.mask;
end
source.dimord = 'pos_freq';

cfg=[];
cfg.parameter = {'stat', 'mask', 'avg'};
source_int = ft_sourceinterpolate(cfg, source, mri);
%%
cmap = flipud(brewermap(64,'RdBu'));
cfgp=[];
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'slice';
cfgp.nslices    =  30;
cfgp.slicerange =  [40 150];
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'stat'; %mask
for k=1:numel(frequency)
    cfgp.frequency = frequency(k);
    ft_sourceplot(cfgp, source_int)
end

% filename = fullfile(['project/3011085.03/', 'analysis', 'source', 'post_cue_pow_stat.m']);
% save(filename, 'source_int', 'source', 'stat', 'frequency','foi');
%% Define ROI's by browsing through Ortho Maps
source = ft_convert_units(source, 'mm');

cfgp=[];
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'ortho';
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'stat';

cfgp.frequency = 50%frequency(1);
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


