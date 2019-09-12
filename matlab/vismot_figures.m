addpath('/project/3011085.03/scripts/RainCloudPlots/tutorial_matlab/')
addpath('/project/3011085.03/scripts/Robust_Statistical_Toolbox/')
subjects = vismot_subjinfo;
alldir = '/project/3011085.03/';
load(fullfile([alldir, 'analysis/source/roi.mat']));
load list

%% FIGURE 2: behavior
% general behavioral measures

for k=1:19
x = load(sprintf('/project/3011085.03/analysis/behaviour/%sbehaviour.mat', list{k}));
performance(k) = sum(x.correct(:))/840;
rttmp = x.rt(x.correct);
rt(k) = mean(rttmp(rttmp~=0));
end
mean(performance)
std(performance)

mean(rt)
std(rt)

for k=1:19
subject= vismot_subjinfo(list{k});
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
ntrl(k) = numel(unique(alldata.data1.trialinfo(:,2)))  + numel(unique(alldata.data2.trialinfo(:,2))) + numel(unique(alldata.data3.trialinfo(:,2))) + numel(unique(alldata.data4.trialinfo(:,2))) + numel(unique(alldata.data5.trialinfo(:,2)));
end
mean(ntrl)
std(ntrl)

% Figure 2a : SIMON EFFECT
filename = [subjects(1).pathname, '/stat_behavior_simon.mat'];
d = load(filename);
effectsize_ms = round(abs(mean(d.avgC-d.avgIC)), 3, 'significant');
effectsize_perc = round(abs(mean(d.avgC./d.avgIC-1)*100),1);
sprintf('behavioral advantage congruent trials: %d ms, %s percent', effectsize_ms, num2str(effectsize_perc))
sprintf('SD = %d, t(%d) = %s, p = %s', round(d.STATS.sd), d.STATS.df, num2str(round(d.STATS.tstat,3, 'significant')), num2str(round(d.P, 3, 'significant')))

figure;
cmap = (brewermap(2,'RdBu'));
f1 = figure;
h2 = raincloud_plot(d.avgIC/1000, 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h1 = raincloud_plot(d.avgC/1000, 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);

legend([h1{1} h2{1}], {'Congruent', 'Incongruent'});
title('reacion times - Simon effect');
set(gca,'Ylim', [-2 5], 'Xlim', [.300 1.100]);
xlabel('Reaction time (s)'); ylabel('Probability density (1/s)')
box off


% Figure 1b: GRATTON EFFECT
cmap = brewermap(2, 'RdBu');
filename = [subjects(1).pathname, '/stat_behavior_gratton.mat'];
d2 = load(filename);
partialeta2 = 0.827;% see /project/3011085.03/analysis/rt/Mats_output_gratton
sprintf('partial eta2 = %s, F(2,17) = %s, p = %s',num2str(round(partialeta2, 3, 'significant')), num2str(round(d2.stat{4,5},2)), num2str(round(d2.stat{4,6},3, 'significant')))


figure;
data{1,1} = d2.avgC_C/1000;
data{2,1} = d2.avgN_C/1000;
data{3,1} = d2.avgIC_C/1000;
data{1,2} = d2.avgC_IC/1000;
data{2,2} = d2.avgN_IC/1000;
data{3,2} = d2.avgIC_IC/1000;
h  = rm_raincloud(data, cmap);
% note that the x and y axes are switched
xlabel('Reaction time (s)')
yticklabels({'previous IC', 'previous N', 'previous C'}) % in reverse order because of switched axes!
xlim([.400 1.100])
legend({'current congruent', 'current incongruent'}, 'Location', 'NorthEast') % how to implement this?



%% FIGURE 3: Post cue power slice plots
alldir = '/project/3011085.03/';
load(fullfile([alldir, 'analysis/source/roi.mat']));
d = load(fullfile([alldir, 'analysis/stat_bf_post.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

effect_size = mean(d.effectsize_largest_cluster); % get average from 'significant clusters' of raw effect (not 1st level T)
effect_size_sd = std(d.effectsize_largest_cluster);
sprintf('p = %s, effect size beta M = %s percent (SD = %s percent), gamma1 M = %s percent (SD = %s percent), gamma2 M = %s percent (SD = %s percent), gamma3 M = %s percent (SD = %s percent).', num2str(round(d.stat.posclusters(1).prob,3 , 'significant')), num2str(round(100*effect_size(2),1)),num2str(round(100*effect_size_sd(2),2)),num2str(round(100*effect_size(3),1)),num2str(round(100*effect_size_sd(3),2)),num2str(round(100*effect_size(4),1)),num2str(round(100*effect_size_sd(4),2)),num2str(round(100*effect_size(5),1)),num2str(round(100*effect_size_sd(5),2)))
d2 = load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']));
effect_size_stratified = mean(d2.effectsize_largest_cluster); % get average from 'significant clusters' of raw effect (not 1st level T)
effect_size_stratified_sd = std(d2.effectsize_largest_cluster);
sprintf('RT-stratified p = %s, effect size alpha M = %s percent (SD = %s percent), effect size beta M = %s percent (SD = %s percent), gamma1 M = %s percent (SD = %s percent), gamma2 M = %s percent (SD = %s percent), gamma3 M = %s percent (SD = %s percent).', num2str(round(d2.stat.posclusters(1).prob,3 , 'significant')), num2str(round(100*effect_size_stratified(1),1)),num2str(round(100*effect_size_stratified_sd(1),2)), num2str(round(100*effect_size_stratified(2),1)),num2str(round(100*effect_size_stratified_sd(2),2)),num2str(round(100*effect_size_stratified(3),1)),num2str(round(100*effect_size_stratified_sd(3),2)),num2str(round(100*effect_size_stratified(4),1)),num2str(round(100*effect_size_stratified_sd(4),2)),num2str(round(100*effect_size_stratified(5),1)),num2str(round(100*effect_size_stratified_sd(5),2)))

% make slice plots of power
cfg=[];
cfg.parameter = {'stat'};
stat_int = ft_sourceinterpolate(cfg, d.stat, mri);

cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
% cfgp.maskparameter = cfgp.funparameter;
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
clim = [5 6 5 8 6];
stat_int.stat(isnan(stat_int.stat)) = 0;
k=1;
for f=stat_int.freq
  cfgp.frequency = f;
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, stat_int);
%   title(sprintf('post-cue power - %d Hz',f))
  k=k+1;
end

f=[];
f.powspctrm = ones(1,10,10);
f.dimord = 'chan_freq_time';
f.time = 1:10;
f.freq = 1:10;
f.label{1} = 'tmp';

cfgp=[];
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.zlim = [-8 8];
figure; ft_singleplotTFR(cfgp, f);




%% FIGURE 4: post cue power ROI

zeroidx = [2 7 14];
meanpower(:,:,1) = d.c_left; meanpower(:,:,2) = d.c_right; meanpower(:,:,3) = d.ic_left; meanpower(:,:,4) = d.ic_right;
stdpower = std(meanpower, [], 3);
meanpower = mean(meanpower,3);
c_ipsi = zeros(19,15);
c_ipsi(:, setdiff(1:15, zeroidx)) = (d.c_left-meanpower)./stdpower;
c_contra = zeros(19,15);
c_contra(:, setdiff(1:15, zeroidx)) = (d.c_right-meanpower)./stdpower;
ic_ipsi = zeros(19,15);
ic_ipsi(:, setdiff(1:15, zeroidx)) = (d.ic_left-meanpower)./stdpower;
ic_contra = zeros(19,15);
ic_contra(:, setdiff(1:15, zeroidx)) = (d.ic_right-meanpower)./stdpower;

figure; 
cnt=1;
for k = setdiff(1:15, zeroidx)
  subplot(3,5,k);
  hold on
  plot([1 2], [ic_ipsi(:,k) c_ipsi(:,k)], '.-', 'color', [0.7 0.7 0.7]);
  plot([1 2], [mean(ic_ipsi(:,k)) mean(c_ipsi(:,k))], 'o-k', 'LineWidth',2)
%   errorbar([1 2], [mean(ic_contra(:,k)) mean(c_contra(:,k))], [std(ic_contra(:,k)) std(c_contra(:,k))])

  plot([3 4], [ic_contra(:,k) c_contra(:,k)], '.-', 'color', [0.7 0.7 0.7]);
  plot([3 4], [mean(ic_contra(:,k)) mean(c_contra(:,k))], 'o-k', 'LineWidth',2)
%   errorbar([3 4], [mean(ic_ipsi(:,k)) mean(c_ipsi(:,k))], [std(ic_ipsi(:,k)) std(c_ipsi(:,k))])
  xlim([0.5 4.5])
  ylim([-2 2])
  cnt=cnt+1;
end



%% Figure 5: pre cue power
alldir = '/project/3011085.03/';
s = load(fullfile([alldir, 'analysis/stat_bf_pre.mat']));
clear meanpower stdpower

zeroidx = [2 7 14];
meanpower(:,:,1) = s.c_left; meanpower(:,:,2) = s.c_right; meanpower(:,:,3) = s.ic_left; meanpower(:,:,4) = s.ic_right;
stdpower = std(meanpower, [], 3);
meanpower = mean(meanpower,3);
c_ipsi = zeros(19,15);
c_ipsi(:, setdiff(1:15, zeroidx)) = (s.c_left-meanpower)./stdpower;
c_contra = zeros(19,15);
c_contra(:, setdiff(1:15, zeroidx)) = (s.c_right-meanpower)./stdpower;
ic_ipsi = zeros(19,15);
ic_ipsi(:, setdiff(1:15, zeroidx)) = (s.ic_left-meanpower)./stdpower;
ic_contra = zeros(19,15);
ic_contra(:, setdiff(1:15, zeroidx)) = (s.ic_right-meanpower)./stdpower;


figure; 
cnt=1;
for k = setdiff(1:15, zeroidx)
  subplot(3,5,k);
  hold on
  plot([1 2], [ic_ipsi(:,k) c_ipsi(:,k)], '.-', 'color', [0.7 0.7 0.7]);
  plot([1 2], [mean(ic_ipsi(:,k)) mean(c_ipsi(:,k))], 'o-k', 'LineWidth',2)
%   errorbar([1 2], [mean(ic_contra(:,k)) mean(c_contra(:,k))], [std(ic_contra(:,k)) std(c_contra(:,k))])
  plot([3 4], [ic_contra(:,k) c_contra(:,k)], '.-', 'color', [0.7 0.7 0.7]);
  plot([3 4], [mean(ic_contra(:,k)) mean(c_contra(:,k))], 'o-k', 'LineWidth',2)
%   errorbar([3 4], [mean(ic_ipsi(:,k)) mean(c_ipsi(:,k))], [std(ic_ipsi(:,k)) std(c_ipsi(:,k))])
  xlim([0.5 4.5])
  
%   title(sprintf('p=%s, p=%s', num2str(round(s.stat.prob(cnt),3)), num2str(round(s.stat.prob(cnt+12),3))))
  cnt=cnt+1;
end

sprintf('After correction there is no statistical significant difference between C and IC trials in the pre-cue window, in any ROI')



%% Figure 6: Coherence
% within-between hemispheres
c = load('/project/3011085.03/analysis/stat_coh_pre.mat');

cmap = brewermap(2,'RdBu');
figure;
subplot(1,2,1)
h1 = raincloud_plot(c.coh_within_C, 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(c.coh_within_IC, 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'Previous congruent', 'Previous incongruent'});
title('Coherence within hemispheres');
set(gca,'Xlim', [0.05 0.3], 'Ylim', [-12 25]);
xlabel('Coherence'); ylabel('Probability (coh^-1)')
box off


subplot(1,2,2)
h3 = raincloud_plot(c.coh_between_C, 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h4 = raincloud_plot(c.coh_between_IC, 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'Previous congruent', 'Previous incongruent'});
title('Coherence between hemispheres');
set(gca,'Xlim', [0.02 0.2], 'Ylim', [-10 20]);
xlabel('Coherence'); ylabel('Probability (coh^-1)')
box off

% individual coherences
stat = c.stat.stat; 
% stat = mean(c.allconnectivity,1); % still not raw coh values-->1st level
% T
coh = zeros(6,6,5);
% alpha
x = ones(6);
x(find(triu(x)))=0;
x(find(x)) = stat(1:c.ncomparisson(1));
coh(:,:,1) = x;

%beta
x = zeros(6);
idx = c.ncomparisson(1)+1;
x(6,3) = stat(idx);
coh(:,:,2) = x;

% gamma1
x = ones(6);
x(find(triu(x)))=0;
idx = sum(c.ncomparisson([1 2]))+1:sum(c.ncomparisson([1 2]))+c.ncomparisson(3);
x(find(x)) = stat(idx);
coh(:,:,3) = x;

% gamma2
x = ones(6);
x(find(triu(x)))=0;
x([3 6],:) = 0; % no motor
x(:,[3 6]) = 0;
idx = sum(c.ncomparisson(1:3))+1:sum(c.ncomparisson(1:3))+c.ncomparisson(4);
x(find(x)) = stat(idx);
coh(:,:,4) = x;

% gamma 3
x = ones(6);
x(find(triu(x)))=0;
idx = sum(c.ncomparisson(1:4))+1:sum(c.ncomparisson(1:4))+c.ncomparisson(5);
x(find(x)) = stat(idx);
coh(:,:,5) = x;

cmap = brewermap(5, 'Accent');
C = [];
C.label = {'left occipital', 'left parietal', 'left motor', 'right occipital', 'right parietal', 'right motor'};
C.channelcmb = ft_channelcombination({'all' 'all'}, C.label);
C.cohspctrm = coh;
C.dimord = 'chan_chan_freq';
C.freq=[10 22 40 60 80];

label = {'L_occ', 'L_par', 'L_mot','R_occ', 'R_par', 'R_mot'};
cmap = brewermap(2, 'RdBu');
addpath('/project/3011085.03/scripts/circularGraph/');
% permute data and labels such that the graph is symmetrical and in the
% right order
order = [2 1 4 5 6 3];
label = label(order);
C.cohspctrm = C.cohspctrm(order, order, :);
for k=1:5
  figure;
  viscircles([0,0], 1, 'color','k'); hold on
  tmpc = C.cohspctrm(:,:,k);
  tmpc(tmpc<0) = 0;
  circularGraph(tmpc, 'Colormap', repmat(cmap(1,:), [6 1]), 'Label', label);
  tmpc = C.cohspctrm(:,:,k);
  tmpc(tmpc>0) = 0;
  circularGraph(abs(tmpc), 'Colormap', repmat(cmap(2,:), [6 1]), 'Label', label);
end


%% Figure S1
d2 = load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']));

cfg=[];
cfg.parameter = {'stat'};
stat_int = ft_sourceinterpolate(cfg, d2.stat, mri);
stat_int.stat(isnan(stat_int.stat)) = 0;

cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
% cfgp.maskparameter = cfgp.funparameter;
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
clim = [5 8 6 8 8];
k=1;
for f=stat_int.freq
  cfgp.frequency = f;
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, stat_int);
%   title(sprintf('post-cue power - %d Hz',f))
  k=k+1;
end

%% Figure S2 ROI plot
s = load(fullfile([alldir, 'analysis/stat_bf_pre.mat'])); l=s.l;r=s.r;

s = d2.stat;
s.stat = d2.stat.stat(:,1);
s.freq = 1;
cmap = brewermap(3, 'Accent');
location = {'L-occipital', 'R-occipital'; 'L-parietal', 'R-parietal'; 'L-motor', 'R-motor'};


cfgp=[];
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'stat';
cfgp.method = 'surface';
cfgp.funcolorlim = [-100, 100];
cfgp.colorbar = 'no';
cfgp.surfdownsample = 20;
cfgp.facecolor = 'brain';

% alpha
ft_sourceplot(cfgp, s);
material dull
p = [1 5 9];
grad=[];
grad.label = {'1'};
yextra =[0 -1.2 0;0 0 0; 0 0 0];
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)] + yextra(k,:);
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 50, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)] + yextra(k,:);
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 50, 'facecolor', [cmap(k,:)])
end
view([0,69]);
% view([-90,0]);camlight

% beta
ft_sourceplot(cfgp, s);
material dull
p = [10];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)]+[0 0 1.5];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(3,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)]+[0 0 1.5];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(3,:)])
end
view([0,69]);
% view([-90,0]);camlight

% gamma1
ft_sourceplot(cfgp, s);
material dull
p = [2 6 11];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  if k==3
    grad.chanpos = grad.chanpos + [0 0 2];
  end
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
    if k==3
    grad.chanpos = grad.chanpos + [0 0 2];
  end
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
% view([-90,0]);camlight

% gamma2
ft_sourceplot(cfgp, s);
material dull
p = [3 7];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  if k==1
    grad.chanpos = grad.chanpos - [0 1 0];
  end
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
  if k==1
    grad.chanpos = grad.chanpos - [0 1.2 0];
  end
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
% view([-90,0]);camlight

% gamma3
ft_sourceplot(cfgp, s);
material dull
p = [4 8 12];
grad=[];
grad.label = {'1'};
yextra1 =[0 -0.5 0;0 0 0; 0 0 0];
yextra2 =[0 -0.5 0;0 0 0; 0 0 0];
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)] + yextra2(k,:);
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)] + yextra1(k,:);
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
% view([-90,0]);camlight

%% S3: pre cue power slice
alldir = '/project/3011085.03/';
s = load(fullfile([alldir, 'analysis/stat_bf_pre.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

source = s.source;
nul=source;
nul.stat = nul.stat*0;

cfgs=[];
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
n=19;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.method = 'analytic';
stat=ft_sourcestatistics(cfgs, source, nul);

cfg=[];
cfg.parameter = {'stat'};
stat_int = ft_sourceinterpolate(cfg, stat, mri);

cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
lim = [5 5 4 5 5];
k=1;
for f=stat_int.freq
cfgp.frequency = f;
cfgp.funcolorlim = [-lim(k) lim(k)];
ft_sourceplot(cfgp, stat_int);
k=k+1;
end


