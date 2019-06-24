addpath('/project/3011085.03/scripts/RainCloudPlots/tutorial_matlab/')
addpath('/project/3011085.03/scripts/Robust_Statistical_Toolbox/')
subjects = vismot_subjinfo;
alldir = '/project/3011085.03/';

%% FIGURE 2 - behavior
filename = [subjects(1).pathname, '/rt/', 'stat_simon.mat'];
d = load(filename);
effectsize_ms = round(abs(mean(d.avgC-d.avgIC)));
effectsize_perc = round(abs(mean(d.avgC./d.avgIC-1)*100),1);
sprintf('behavioral advantage congruent trials: %d ms, %s percent', effectsize_ms, num2str(effectsize_perc))
sprintf('SD = %d, t(%d) = %s, p = %s', round(d.STATS.sd), d.STATS.df, num2str(round(d.STATS.tstat,3, 'significant')), num2str(round(d.P, 3, 'significant')))

% Figure 1a : SIMON EFFECT
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
filename = [subjects(1).pathname, '/rt/', 'stat_gratton.mat'];
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



%% FIGURE 3 - Post cue power slice plots
alldir = '/project/3011085.03/';
load(fullfile([alldir, 'analysis/source/roi.mat']));
d = load(fullfile([alldir, 'analysis/stat_bf_post.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

effect_size = mean(d.effectsize_largest_cluster); % get average from 'significant clusters' of raw effect (not 1st level T)
effect_size_sd = std(d.effectsize_largest_cluster);
sprintf('effect size M = %s percent (SD = %s percent), p = %d', num2str(round(100*effect_size,1)),num2str(round(100*effect_size_sd,1)), d.stat.posclusters(1).prob)

d2 = load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']));
effect_size_stratified = mean(d2.effectsize_largest_cluster); % get average from 'significant clusters' of raw effect (not 1st level T)
effect_size_stratified_sd = std(d2.effectsize_largest_cluster);
sprintf('RT-stratified effect size M = %s percent (SD = %s percent), p = %s', num2str(round(100*effect_size_stratified,1)),num2str(round(100*effect_size_stratified_sd,1)), num2str(round(d2.stat.posclusters(1).prob, 3, 'significant')))

% make slice plots of power
cfg=[];
cfg.parameter = 'stat';
stat_int = ft_sourceinterpolate(cfg, d2.stat, mri);

cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.maskparameter = cfgp.funparameter;
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
clim = [5 6 5 8 6];
k=1;
for f=stat_int.freq
  cfgp.frequency = f;
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, stat_int);
  title(sprintf('post-cue power - %d Hz',f))
  k=k+1;
end


%% FIGURE 4 - post cue power ROI + bar plot
% Make bar graph of defined of power in defined ROIs and FOIs
% find indices of roi locations

idx_left = d.l;
idx_right = d.r;

cmap = brewermap(5, 'Accent');

figure;
% occipital
subplot(1,3,1);
occ = [mean(d.effectsize_roi_left(:,1),1), mean(d.effectsize_roi_left(:,2),1) mean(d.effectsize_roi_left(:,3),1), mean(d.effectsize_roi_left(:,4),1); mean(d.effectsize_roi_right(:,1),1), mean(d.effectsize_roi_right(:,2),1) mean(d.effectsize_roi_right(:,3),1), mean(d.effectsize_roi_right(:,4),1)];
% occ = [d.stat.stat(idx_left(1),1), d.stat.stat(idx_left(2),3) d.stat.stat(idx_left(3),4), d.stat.stat(idx_left(4),5); d.stat.stat(idx_right(1),1), d.stat.stat(idx_right(2),3) d.stat.stat(idx_right(3),4), d.stat.stat(idx_right(4),5)];
b = bar(occ);
tmpcmap = cmap([1 3:5],:);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1', 'gamma2', 'gamma3'}, 'location', 'northwest')
title('occipital')

% parietal
subplot(1,3,2);
par = [mean(d.effectsize_roi_left(:,5),1), mean(d.effectsize_roi_left(:,6),1) mean(d.effectsize_roi_left(:,7),1), mean(d.effectsize_roi_left(:,8),1); mean(d.effectsize_roi_right(:,5),1), mean(d.effectsize_roi_right(:,6),1) mean(d.effectsize_roi_right(:,7),1), mean(d.effectsize_roi_right(:,8),1)];
% par = [stat.stat(idx_left(5),1), stat.stat(idx_left(6),3) stat.stat(idx_left(7),4), stat.stat(idx_left(8),5); stat.stat(idx_right(5),1), stat.stat(idx_right(6),3) stat.stat(idx_right(7),4), stat.stat(idx_right(8),5)];
b = bar(par);
tmpcmap = cmap([1 3:5],:);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1', 'gamma2', 'gamma3'}, 'location', 'northwest')
title('parietal')

% motor
subplot(1,3,3);
% mot = [stat.stat(idx_left(9),1), stat.stat(idx_left(10),2) stat.stat(idx_left(11),3), stat.stat(idx_left(12),5); stat.stat(idx_right(9),1), stat.stat(idx_right(10),2) stat.stat(idx_right(11),3), stat.stat(idx_right(12),5)];
mot = [mean(d.effectsize_roi_left(:,9),1), mean(d.effectsize_roi_left(:,10),1) mean(d.effectsize_roi_left(:,11),1), mean(d.effectsize_roi_left(:,12),1); mean(d.effectsize_roi_right(:,9),1), mean(d.effectsize_roi_right(:,10),1) mean(d.effectsize_roi_right(:,11),1), mean(d.effectsize_roi_right(:,12),1)];

b = bar(mot);
tmpcmap = cmap([1:3 5],:);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'beta', 'gamma1', 'gamma3'}, 'location', 'northeast')
title('motor')

%% or per frequency
% Make bar graph of defined of power in defined ROIs and FOIs
% find indices of roi locations
for k=1:size(roi,1)-1
    idx_left(k) = find_dipoleindex(d2.stat, roi{k+1,3});
    idx_right(k) = find_dipoleindex(d2.stat, roi{k+1,4});
end


cmap = brewermap(3, 'Accent');
freqs = {'8-12', '14-30', '30-50', '50-70', '70-90'};

figure;
for f=1:5
subplot(1,5,f);
a = [stat.stat(idx_left(1),f), stat.stat(idx_left(5),f) stat.stat(idx_left(9),f); stat.stat(idx_right(1),f), stat.stat(idx_right(5),f) stat.stat(idx_right(9),f)];
b = bar(a);
for k=1:3
  set(b(k),'FaceColor', cmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
if f==1; legend({'occipital', 'parietal', 'motor'}), end
title(sprintf('%s Hz', freqs{f}))
end

% ROI plot
s = stat;
s.stat = s.stat(:,1);
s.freq = 1;
load('/project/3011085.03/analysis/stat_bf_pre.mat', 'l', 'r');
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
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 50, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 50, 'facecolor', [cmap(k,:)])
end
view([0,69]);
view([-90,0]);camlight

% beta
ft_sourceplot(cfgp, s);
material dull
p = [10];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)]+[0 0 2];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(3,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)]+[0 0 2];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(3,:)])
end
view([0,69]);
view([-90,0]);camlight

% gamma1
ft_sourceplot(cfgp, s);
material dull
p = [2 6 11];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
view([-90,0]);camlight

% gamma2
ft_sourceplot(cfgp, s);
material dull
p = [3 7];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
view([-90,0]);camlight

% gamma3
ft_sourceplot(cfgp, s);
material dull
p = [4 8 12];
grad=[];
grad.label = {'1'};
for k=1:length(p)
  grad.chanpos=[s.pos(l(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(k,:)])
end
view([0,69]);
view([-90,0]);camlight

%% pre cue power
alldir = '/project/3011085.03/';
s = load(fullfile([alldir, 'analysis/stat_bf_pre.mat']));
% OPTION 1
cmap = brewermap(2, 'RdBu');
figure;
for k=1:12
d{k,1} = s.lpow(:,k);
d{k,2} = s.rpow(:,k);
end
h = rm_raincloud(d, cmap);
for k=1:numel(h.l)
  h.l(k).Color = 'none';
end
for k=1:numel(h.m)
    h.m(k).MarkerFaceColor = 'none';
  h.m(k).MarkerEdgeColor = 'none';
end

% OPTION 2
figure;
cmap2 = brewermap(2, 'RdBu');
for k=1:size(s.lpow,2)
  subplot(3,4,k)
  h1{k} = raincloud_plot(s.lpow(:,k), 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
  h2{k} = raincloud_plot(s.rpow(:,k), 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
end

% OPTION 3
  diff = s.lpow-s.rpow;
  addpath /project/3011085.03/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(diff', 'showScatter', true, 'scatterMarker', '.'); 
load(fullfile([alldir, 'analysis/source/roi.mat']));
cmap = brewermap(5, 'Accent');
for k=1:12
%   c{k} = cmap(s.freq_idx(k),:);
x{1}(k,:) = [cmap(s.freq_idx(k),:)];
end

% OPTION 4
meanl = mean(s.lpow,1);
meanr = mean(s.rpow,1);
figure;
tmpcmap = cmap([1 3:5],:);
subplot(1,3,1)
b = bar([s.lpow(1:4); s.rpow(1:4)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
subplot(1,3,2)
b = bar([s.lpow(5:8); s.rpow(5:8)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
subplot(1,3,3)
tmpcmap = cmap([1:3 5],:);
b = bar([s.lpow(9:12); s.rpow(9:12)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});


%% Coherence
% within-between hemispheres
c = load('/project/3011085.03/analysis/stat_coh_pre.mat')

cmap = flipud(brewermap(2,'RdBu'));
f1 = figure;
h1 = raincloud_plot(c.zx_within, 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(c.zx_between, 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'Intrahemispheric connectivity', 'Interhemispheric connectivity'});
title('Coherence within vs between hemispheres');
set(gca,'Xlim', [-0.03 0.03], 'Ylim', [-30 60]);
xlabel('Coherence congruent-incongruent (T)'); ylabel('Probability (a.u.)')
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

% FT way
% make a parameter containing width of the arrows (in pixels) according to
% the coh value
C.width = round(2*C.cohspctrm);
h = figure;
viscircles([0,0], 1, 'color','k')
load('vismot_parcels', 'layout3', 'plotlayoptions3');
layout3.label = C.label;
layout3.pos = [layout3.pos(4:6,:); layout3.pos(1:3,:)];
% plotlayoptions3{10} = [plotlayoptions3{10}(4:6); plotlayoptions3{10}(1:3)];
plotlayoptions3{8}=zeros(1,6);

cfgp=[];
cfgp.layout=layout3;
cfgp.widthparam = 'width';
cfgp.newfigure = 'no';
cfgp.colormap = flipud(brewermap(2, 'RdBu'));
for k=1:numel(C.freq)
  figure; hold on; 
%   viscircles([0,0], 1, 'color','k')
  cfgp.colormap = cmap(k,:);
  cfgp.foi = C.freq(k);
  ft_topoplotCC(cfgp, C);
  ft_plot_layout(layout3, plotlayoptions3{:});
end
% ft_plot_layout(layout3, plotlayoptions3{:});


% Circular Graph
label = {'L_occ', 'L_par', 'L_mot','R_occ', 'R_par', 'R_mot'};
cmap = brewermap(5, 'Accent');
addpath('/project/3011085.03/scripts/circularGraph/');
for k=1:5
  figure;
  viscircles([0,0], 1, 'color','k'); hold on
  circularGraph(C.cohspctrm(:,:,k), 'Colormap', repmat(cmap(k,:), [6 1]), 'Label', label);
end

% or 
label = {'L_occ', 'L_par', 'L_mot','R_occ', 'R_par', 'R_mot'};
cmap = brewermap(2, 'RdBu');
addpath('/project/3011085.03/scripts/circularGraph/');
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




% legend({'alpha', 'beta', 'gamma1', 'gamma2', 'gamma3'}) % doesnt work