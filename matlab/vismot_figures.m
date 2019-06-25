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
% raincloud graph of ipsi vs contra ROI pow diff (C./IC)
idx_left = d2.l;
idx_right = d2.r;

pow_ipsi = zeros(19,15);
pow_contra = zeros(19,15);
zeroidx = [2 7 14];
pow_ipsi(:, setdiff(1:15, zeroidx)) = d2.effectsize_roi_left;
pow_contra(:, setdiff(1:15, zeroidx))  = d2.effectsize_roi_right;

cmap = brewermap(2, 'RdBu');
figure;
for k = setdiff(1:15, zeroidx)
    subplot(3,5,k)
  h1{k} = raincloud_plot(pow_ipsi(:,k), 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
  h2{k} = raincloud_plot(pow_contra(:,k), 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
f1 = ksdensity(pow_ipsi(:,k));
f2 = ksdensity(pow_contra(:,k));
   maxy = max([f1 f2]);
   set(gca, 'Ylim', [-maxy maxy])
end

% ROI plot
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
  grad.chanpos=[s.pos(l(p(k)),:)]+[0 0 1.5];
  grad.elecpos=grad.chanpos;
  ft_plot_sens(grad,'elecshape','point', 'elecsize', 40, 'facecolor', [cmap(3,:)])
  
  grad.chanpos=[s.pos(r(p(k)),:)]+[0 0 1.5];
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
view([-90,0]);camlight

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

pow_ipsi = zeros(19,15);
pow_contra = zeros(19,15);
zeroidx = [2 7 14];
pow_ipsi(:, setdiff(1:15, zeroidx)) = s.effectsize_roi_left;
pow_contra(:, setdiff(1:15, zeroidx))  = s.effectsize_roi_right;

figure;
cmap = brewermap(2, 'RdBu');
for k=setdiff(1:15, zeroidx)
  subplot(3,5,k)
  h1{k} = raincloud_plot(pow_ipsi(:,k), 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
  h2{k} = raincloud_plot(pow_contra(:,k), 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
   f1 = ksdensity(pow_ipsi(:,k));
f2 = ksdensity(pow_contra(:,k));
   maxy = max([f1 f2]);
   set(gca, 'Ylim', [-maxy maxy])
end



%% Coherence
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
