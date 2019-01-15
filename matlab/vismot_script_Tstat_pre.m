load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = 'statResp';
frequency = [4:2:30 34:4:94];
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%ssource3d4mm_pre_%03d.mat', list{m}, k)]);
    dum = load(d, whichstat);
    tmp(m) = dum.(whichstat);
  end
  dat(:, cnt+1,1:n) = cat(2,tmp.stat);
  clear tmp dum
  cnt=cnt+1;
end
foi = { 'alpha', 10, 10
        'beta2', [24 26], 25
        'gamma1', [58 62], 60
        'gamma2', [78 82], 80};
    
dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = [10 25 60 80];
source.stat(:,:,1) = dat(:,:,find(frequency==10));
source.stat(:,:,2) = nanmean(dat(:,:,find(frequency==24):find(frequency==26)),3);
source.stat(:,:,3) = nanmean(dat(:,:,find(frequency==58):find(frequency==62)),3);
source.stat(:,:,4) = nanmean(dat(:,:,find(frequency==78):find(frequency==82)),3);

source.avg = squeeze(nanmean(source.stat,1));

nul=source;
nul.stat(:)=0;
%% look at FOIs, determined in post-cue window
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
for k=1:4
    source.stat(:,k) = stat{k}.stat;
    source.mask(:,k) = stat{k}.mask;
end
source.dimord = 'pos_freq';

cfg=[];
cfg.parameter = {'stat', 'mask', 'avg'};
source_int = ft_sourceinterpolate(cfg, source, mri);

cmap = flipud(brewermap(64,'RdBu'));
cfgp=[];
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'slice';
cfgp.nslices    =  30;
cfgp.slicerange =  [40 150];
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'stat'; %mask
for k=1:numel(source.freq)
    cfgp.frequency = source.freq(k);
    ft_sourceplot(cfgp, source_int)
end



%% Look at spectra in ROIs
s = rmfield(source, {'avg', 'stat'});
s.freq = frequency;
s.stat = dat;
s.avg = squeeze(nanmean(s.stat,1));
s.dimord = 'rpt_pos_freq';

% load ROIs
load('/project/3011085.03/analysis/source/roi.mat', 'ROI');
for k=1:3
    idx_left(k) = find_dipoleindex(source, ROI{k,2});
    idx_right(k) = find_dipoleindex(source, ROI{k,3});
end

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'no';
cfgs.numrandomization = 1000;
cfgs.clusteralpha = 0.05;
cfgs.correcttail = 'prob';

s.inside(:)=0;
s.inside([idx_left, idx_right]) = 1;
nul=s;
nul.stat(:)=0;
s = ft_sourcestatistics(cfgs, s, nul);

% plot 2nd level tstat against zero
% for FOIs
figure;
% occipital
k=1;
subplot(1,3,k);
y = [nanmean(s.stat(idx_left(k), find(s.freq==10))), nanmean(s.stat(idx_left(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_left(k), find(s.freq==42))), nanmean(s.stat(idx_left(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_left(k),[find(s.freq==78) find(s.freq==82)])); nanmean(s.stat(idx_right(k), find(s.freq==10))), nanmean(s.stat(idx_right(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_right(k), find(s.freq==42))), nanmean(s.stat(idx_right(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_right(k),[find(s.freq==78) find(s.freq==82)]))]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
% legend({'alpha', 'beta 20-30', 'gamma 42', 'gamma 50-70', 'gamma 70-90'})
title('occipital')

% parietal
k=2;
subplot(1,3,k);
y = [nanmean(s.stat(idx_left(k), find(s.freq==10))), nanmean(s.stat(idx_left(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_left(k), find(s.freq==42))), nanmean(s.stat(idx_left(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_left(k),[find(s.freq==78) find(s.freq==82)])); nanmean(s.stat(idx_right(k), find(s.freq==10))), nanmean(s.stat(idx_right(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_right(k), find(s.freq==42))), nanmean(s.stat(idx_right(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_right(k),[find(s.freq==78) find(s.freq==82)]))]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'beta 20-30', 'gamma 42', 'gamma 50-70', 'gamma 70-90'}, 'location', 'northwest')
title('parietal')

% motor
k=3;
subplot(1,3,k);
y = [nanmean(s.stat(idx_left(k), find(s.freq==10))), nanmean(s.stat(idx_left(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_left(k), find(s.freq==42))), nanmean(s.stat(idx_left(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_left(k),[find(s.freq==78) find(s.freq==82)])); nanmean(s.stat(idx_right(k), find(s.freq==10))), nanmean(s.stat(idx_right(k),[find(s.freq==24) find(s.freq==26)])), nanmean(s.stat(idx_right(k), find(s.freq==42))), nanmean(s.stat(idx_right(k),[find(s.freq==58) find(s.freq==62)])), nanmean(s.stat(idx_right(k),[find(s.freq==78) find(s.freq==82)]))]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
% legend({'alpha', 'beta 20-30', 'gamma 42', 'gamma 50-70', 'gamma 70-90'})
title('motor')


figure;
subplot(1,3,1);
plot(s.freq, s.stat(idx_left(1),:)); hold on; plot(s.freq, s.stat(idx_right(1),:));
title('occipital')
legend('left', 'right');
subplot(1,3,2);
plot(s.freq, s.stat(idx_left(2),:)); hold on; plot(s.freq, s.stat(idx_right(2),:));
title('parietal')
legend('left', 'right');
subplot(1,3,3);
plot(s.freq, s.stat(idx_left(3),:)); hold on; plot(s.freq, s.stat(idx_right(3),:));
title('motor')
legend('left', 'right');
suptitle('T-values in ROIs (specified in post-cue window)')
