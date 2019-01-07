preorpost = input('pre or post-cue?');

load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
whichstat = 'statResp';
frequency = [4:2:30 34:4:94];%[6 10 16 24 26 58 62 78 82];%
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 0;
for k = frequency
for m = 1:n
d = fullfile([datadir, sprintf('%ssource3d4mm_%s_%03d.mat', list{m}, preorpost, k)]);
dum = load(d, whichstat);
tmp(m) = dum.(whichstat);
end
dat(:, cnt+1,1:n) = cat(2,tmp.stat);
clear tmp dum
cnt=cnt+1;
end
dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.stat=dat;
source.freq=frequency;


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
cfgs.correctm = 'cluster';
cfgs.numrandomization = 1000;
cfgs.clusteralpha = 0.01;
cfgs.correcttail = 'no';
source.inside(:)=0;
source.inside([idx_left, idx_right]) = 1;
nul=source;
nul.stat(:)=0;
stat = ft_sourcestatistics(cfgs, source, nul);

stat_left = stat.stat(idx_left,:);
stat_right = stat.stat(idx_right,:);
% figure; 
subplot(1,3,3);
plot(frequency, stat_left(1,:),'b'); %hold on; plot(frequency, stat_right(1,:));
title('occipital'); xlabel('frequency'); ylabel('T-value');ylim([-8 8])
% subplot(1,3,2)
% plot(frequency, stat_left(2,:)); hold on; plot(frequency, stat_right(2,:));
% title('parietal'); xlabel('frequency');
% subplot(1,3,3);
% plot(frequency, stat_left(3,:)); hold on; plot(frequency, stat_right(3,:));
% title('motor'); xlabel('frequency'); ylabel('T-value');
% subplot(1,3,2); legend({'left', 'right'}, 'location', 'south')