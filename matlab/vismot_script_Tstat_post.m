load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = 'statResp';
frequency = [6 10 16 24 26 58 62 78 82];%[10 16 24 26 48:4:92];%[4:2:30 40:4:100];
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
        'gamma1', [58 62], 60
        'gamma2', [78 82], 80};
frequency = [6 10 16 25 60 80];
    
dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = [6 10 16 25 60 80];
source.stat(:,:,1:3) = dat(:,:,1:3);
source.stat(:,:,4) = nanmean(dat(:,:,4:5),3);
source.stat(:,:,5) = nanmean(dat(:,:,6:7),3);
source.stat(:,:,6) = nanmean(dat(:,:,8:9),3);
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
cfgs.correcttail = 'prob';
for k=1:size(foi,1)
    cfgs.frequency = source.freq(k);
    stat{k} = ft_sourcestatistics(cfgs, source, nul);
end

source = rmfield(source, 'stat');
for k=1:6
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
cfgp.maskparameter = 'mask'; %mask
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

cfgp.frequency = frequency(1);
cfgp.location='max';
ft_sourceplot(cfgp, source);


