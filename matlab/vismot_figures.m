subjects = vismot_subjinfo;
alldir = '/project/3011085.03/';
addpath([alldir 'scripts/RainCloudPlots/tutorial_matlab/'])
addpath([alldir 'scripts/Robust_Statistical_Toolbox/'])
addpath([alldir 'scripts/circularGraph/']);
load(fullfile([alldir, 'analysis/roi.mat']));
load list

%% FIGURE 2: behavior
% general behavioral measures

for k=1:19
    x = load([alldir sprintf('analysis/behaviour/%sbehaviour.mat', list{k})]);
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
partialeta2 = 0.827;% see [alldir 'analysis/rt/Mats_output_gratton']
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
d = load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

% make slice plots of power
cfg=[];
cfg.parameter = {'stat', 'mask'};
for k=1:6
    stat_int{k} = ft_sourceinterpolate(cfg, d.stat{k}, mri);
end


cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.maskparameter = 'mask2';
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
cfgp.opacitylim = [-0.1 0.1];
clim = [5 6 8 6 8 6];
tmp = {'theta', 'alpha', 'beta', 'gamma1', 'gamma2', 'gamma3'};

for k=1:6
  s=stat_int{k};
  s.mask2 = s.stat;
  xmin = min(s.stat); xmax = max(s.stat);
  s.mask2(s.mask2<0.3*xmax & s.mask2>0.3*xmin & s.inside(:)==1)=0;
  s.mask2(s.mask2>0|s.mask<0)=1;
  s.mask2(isnan(s.mask2))=0;
  ix=find(s.inside(:)==0 & s.anatomy(:)~=0); s.mask2(ix) = 0; % this line makes sure the anatomy doesn't disappear
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, s);
  saveas(gcf, sprintf('3_postpow_slice_%s_str', tmp{k}), 'png')
end


f=[];
f.powspctrm = ones(1,10,10);
f.dimord = 'chan_freq_time';
f.time = 1:10;
f.freq = 1:10;
f.label{1} = 'tmp';

Uclim = unique(clim);
figure;
for k=1:numel(Uclim)
  cfgp=[];
  cfgp.colormap = flipud(brewermap(64, 'RdBu'));
  cfgp.zlim = [-Uclim(k) Uclim(k)];
  ft_singleplotTFR(cfgp, f);
  saveas(gcf, sprintf('colorbar_%d',Uclim(k)), 'eps');
end


%% Figure 4: pre cue power whole brain
d = load(fullfile([alldir, 'analysis/stat_bf_pre_wholebrain.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

% make slice plots of power
cfg=[];
cfg.parameter = {'stat', 'mask'};
for k=1:6
    stat_int{k} = ft_sourceinterpolate(cfg, d.stat{k}, mri);
end

tmp = {'theta', 'alpha', 'beta', 'gamma1', 'gamma2', 'gamma3'};
cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.maskparameter = 'mask2';
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
cfgp.opacitylim = [-0.1 0.1];
clim = [5 5 5 5 5 5];
for k=1:6
  s=stat_int{k};
  s.mask2 = s.stat;
  xmin = min(s.stat); xmax = max(s.stat);
  s.mask2(s.mask2<0.3*xmax & s.mask2>0.3*xmin & s.inside(:)==1)=0;
  s.mask2(s.mask2>0|s.mask<0)=1;
  s.mask2(isnan(s.mask2))=0;
  ix=find(s.inside(:)==0 & s.anatomy(:)~=0); s.mask2(ix) = 0; % this line makes sure the anatomy doesn't disappear
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, s);
  saveas(gcf, sprintf('4_prepow_slice_%s', tmp{k}), 'png')
end


f=[];
f.powspctrm = ones(1,10,10);
f.dimord = 'chan_freq_time';
f.time = 1:10;
f.freq = 1:10;
f.label{1} = 'tmp';

Uclim = unique(clim);
figure;
for k=1:numel(Uclim)
  cfgp=[];
  cfgp.colormap = flipud(brewermap(64, 'RdBu'));
  cfgp.zlim = [-Uclim(k) Uclim(k)];
  ft_singleplotTFR(cfgp, f);
  saveas(gcf, sprintf('colorbar_%d',Uclim(k)), 'eps');
end


%% Figure 5: Coherence 
%%%%%%%%%%%%%%%%%
% RAW COHERENCE %
%%%%%%%%%%%%%%%%%
c = load([alldir 'analysis/stat_coh_pre.mat']);
freqs = {'theta', 'alpha', 'beta', 'gamma2', 'gamma3'};
iscoh=zeros(numel(freqs),1);
for k=1:numel(freqs), cohmean{k} = zeros(7,7); cohstd{k} = zeros(7,7); end
for k=1:numel(c.idx_sign_uncorrected)
    idx = c.idx_sign_uncorrected(k);
    e=c.effect(idx);
    
    % find ROI
    if strfind(e.roi1, 'occipital'), roi1 = 1;
    elseif strfind(e.roi1, 'parietal'), roi1 = 2;
    elseif strfind(e.roi1, 'motor'), roi1 = 3; end
    if strfind(e.roi2, 'occipital'), roi2 = 1;
    elseif strfind(e.roi2, 'parietal'), roi2 = 2;
    elseif strfind(e.roi2, 'motor'), roi2 = 3; end

    % find lateralization
    if strfind(e.roi1, 'ipsi'), l1 = 0;
    elseif strfind(e.roi1, 'contra'), l1 = 3; end
    if strfind(e.roi2, 'ipsi'), l2 = 0;
    elseif strfind(e.roi2, 'contra'), l2 = 3; end
    
    for ii=1:numel(freqs), if strcmp(freqs{ii}, e.freq), fidx = ii; end, end
    
    cohmean{fidx}(l1+roi1, l2+roi2) = 100*mean(c.allcohC(:,idx)./c.allcohIC(:,idx)-1);
    cohstd{fidx}(l1+roi1, l2+roi2) = 100*std(c.allcohC(:,idx)./c.allcohIC(:,idx)-1);
    iscoh(fidx)=1;
end
iscoh = find(iscoh);

warning('please enable line 74 in circularGraph.m')
cmap = brewermap(2, 'RdBu');
transparency = 0.3;
cmap2=[cmap [transparency; transparency]];
% permute data and labels such that the graph is symmetrical and in the
label = {'occ_ipsi', 'par_ipsi', 'mot_ipsi','occ_contra', 'par_contra', 'mot_contra', 'fron'};
order = [2 1 4 5 6 7 3];
label=label(order);

maxnumber = abs(cat(3,cohmean{:}));
maxnumber = max(maxnumber(:));
maxnumberstd = (cat(3,cohstd{:}));
maxnumberstd = max(maxnumberstd(:));
for k=1:numel(iscoh)
    figure;
    viscircles([0,0], 1, 'color','k'); hold on
    suptitle(sprintf('mean %s, sign uncor',freqs{iscoh(k)}))
    vismot_circularGraph(cohmean{iscoh(k)}(order, order),'std', cohstd{iscoh(k)}(order, order), 'Colormap', repmat(cat(3,cmap(1,:), cmap(2,:)), [7 1 1]), 'Label', label, 'maxnumber', maxnumber, 'maxnumberstd', maxnumberstd);
end


%%%%%%%%%%%%
% T-VALUES %
%%%%%%%%%%%%
iscoh=zeros(numel(freqs),1);
for k=1:numel(freqs), coh{k} = zeros(7,7); end
for k=1:numel(c.idx_sign_uncorrected)
    idx = c.idx_sign_uncorrected(k);
    e=c.effect(idx);
    
    % find ROI
    if strfind(e.roi1, 'occipital'), roi1 = 1;
    elseif strfind(e.roi1, 'parietal'), roi1 = 2;
    elseif strfind(e.roi1, 'motor'), roi1 = 3; end
    if strfind(e.roi2, 'occipital'), roi2 = 1;
    elseif strfind(e.roi2, 'parietal'), roi2 = 2;
    elseif strfind(e.roi2, 'motor'), roi2 = 3; end

    % find lateralization
    if strfind(e.roi1, 'ipsi'), l1 = 0;
    elseif strfind(e.roi1, 'contra'), l1 = 3; end
    if strfind(e.roi2, 'ipsi'), l2 = 0;
    elseif strfind(e.roi2, 'contra'), l2 = 3; end
    
    for ii=1:numel(freqs), if strcmp(freqs{ii}, e.freq), fidx = ii; end, end
    
    coh{fidx}(l1+roi1, l2+roi2) = e.stat;
    iscoh(fidx)=1;
end
iscoh = find(iscoh);

warning('please enable line 74 in circularGraph.m')
cmap = brewermap(2, 'RdBu');
% permute data and labels such that the graph is symmetrical and in the
label = {'occ_ipsi', 'par_ipsi', 'mot_ipsi','occ_contra', 'par_contra', 'mot_contra', 'fron'};
order = [2 1 4 5 6 7 3];
label=label(order);
maxnumber = abs(cat(3,coh{:}));
maxnumber = max(maxnumber(:));
for k=1:numel(iscoh)
    figure;
    viscircles([0,0], 1, 'color','k'); hold on
    suptitle(sprintf('%s, sign uncor',freqs{iscoh(k)}))
    vismot_circularGraph(coh{iscoh(k)}(order, order), 'Colormap', repmat(cat(3,cmap(1,:), cmap(2,:)), [7 1 1]), 'Label', label,'maxnumber', maxnumber);
end

%% Figure S1: Post cue power (not stratified)
d = load(fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']));
mri = ft_read_mri('single_subj_T1_1mm.nii');

% make slice plots of power
cfg=[];
cfg.parameter = {'stat', 'mask'};
for k=1:6
    stat_int{k} = ft_sourceinterpolate(cfg, d.stat{k}, mri);
end


cmap = flipud(brewermap(64, 'RdBu'));
cfgp = [];
cfgp.funparameter = 'stat';
cfgp.method = 'slice';
cfgp.funcolormap = cmap;
cfgp.maskparameter = 'mask2';
cfgp.slicerange = [70 150];
cfgp.nslices = 16;
cfgp.opacitylim = [-0.1 0.1];
clim = [5 6 8 6 8 6];
tmp = {'theta', 'alpha', 'beta', 'gamma1', 'gamma2', 'gamma3'};

for k=1:6
  s=stat_int{k};
  s.mask2 = s.stat;
  xmin = min(s.stat); xmax = max(s.stat);
  s.mask2(s.mask2<0.3*xmax & s.mask2>0.3*xmin & s.inside(:)==1)=0;
  s.mask2(s.mask2>0|s.mask<0)=1;
  s.mask2(isnan(s.mask2))=0;
  ix=find(s.inside(:)==0 & s.anatomy(:)~=0); s.mask2(ix) = 0; % this line makes sure the anatomy doesn't disappear
  cfgp.funcolorlim = [-clim(k) clim(k)];
  ft_sourceplot(cfgp, s);
  saveas(gcf, sprintf('3_postpow_slice_%s_str', tmp{k}), 'png')
end


f=[];
f.powspctrm = ones(1,10,10);
f.dimord = 'chan_freq_time';
f.time = 1:10;
f.freq = 1:10;
f.label{1} = 'tmp';

Uclim = unique(clim);
figure;
for k=1:numel(Uclim)
  cfgp=[];
  cfgp.colormap = flipud(brewermap(64, 'RdBu'));
  cfgp.zlim = [-Uclim(k) Uclim(k)];
  ft_singleplotTFR(cfgp, f);
  saveas(gcf, sprintf('colorbar_%d',Uclim(k)), 'eps');
end


%% Figure S2: pre cue power ROI
d = load(fullfile([alldir, 'analysis/stat_bf_pre.mat']));
clear meanpower stdpower c ic

c = (d.c-squeeze(mean(d.raw)))./squeeze(std(d.raw));
ic = (d.ic-squeeze(mean(d.raw)))./squeeze(std(d.raw));
Y = [-0.5 1; -0.5 1; -0.5 1; -1 2; -1 2; -1 2; -1 2;-1 2; -2 4; -2 4; -1 2; -1 2;-1 2; -2 4; -1 2; -1 2; -1 2];
X = [2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1 1.5 1.5 1.5 1.5 1 1 1.5 1.5];
figure
for k=1:17
    subplot(3,6,k+1)
cmap = (brewermap(2,'RdBu'));
h2 = raincloud_plot(c(:,k), 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h1 = raincloud_plot(ic(:,k), 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
    'box_col_match', 0);
% legend([h1{1} h2{1}], {'Congruent', 'Incongruent'});
% title('reacion times - Simon effect');
set(gca,'Ylim', Y(k,:), 'Xlim', [-X(k) X(k)]);
% xlabel('raw power (a.u.)');
ylabel('Probability density (1/s)')
box off
end




stdpower = squeeze(std(d.raw, [], 1));
meanpower = squeeze(mean(d.raw,1));
c=(d.c-meanpower)./stdpower;
ic=(d.ic-meanpower)./stdpower;
figure;
for k = 1:17
hold on
plot([2*(k-1)+1 2*k], [ic(:,k) c(:,k)], '.-', 'color', [0.7 0.7 0.7]);
plot([2*(k-1)+1 2*k], [mean(ic(:,k)) mean(c(:,k))], 'o-k', 'LineWidth',2)
end

% or

c = (d2.c(:,1:17)+d2.c(:,18:end))./2;
ic = (d2.ic(:,1:17)+d2.ic(:,18:end))./2;
figure; hold on
m=repmat(max([c;ic]),19, 1);
boxplot((c-ic)./m)

sprintf('After correction there is no statistical significant difference between C and IC trials in the pre-cue window, in any ROI')

%% Figure S3: coherence of connections not significant (uncorrected)
c = load([alldir 'analysis/stat_coh_pre.mat']);
freqs = {'theta', 'alpha', 'beta', 'gamma2', 'gamma3'};

%%%%%%%%%%%%%%%%%
% RAW COHERENCE %
%%%%%%%%%%%%%%%%%
for k=1:numel(freqs), cohmean2{k} = zeros(7,7);cohstd2{k} = zeros(7,7); end
ncomp2=zeros(7,7,5);
idx_nonsig = setdiff(1:numel(c.effect), c.idx_sign_uncorrected');
for k=1:numel(idx_nonsig)
    idx=idx_nonsig(k);
    e=c.effect(idx);
    
    % find ROI
    if strfind(e.roi1, 'occipital'), roi1 = 1;
    elseif strfind(e.roi1, 'parietal'), roi1 = 2;
    elseif strfind(e.roi1, 'motor'), roi1 = 3; 
    elseif strfind(e.roi1, 'frontal'), roi1 = 4; end
    if strfind(e.roi2, 'occipital'), roi2 = 1;
    elseif strfind(e.roi2, 'parietal'), roi2 = 2;
    elseif strfind(e.roi2, 'motor'), roi2 = 3; 
    elseif strfind(e.roi2, 'frontal'), roi2 = 4; end

    % find lateralization
    if strfind(e.roi1, 'ipsi'), l1 = 0;
    elseif strfind(e.roi1, 'contra'), l1 = 3; 
    elseif strfind(e.roi1, 'midline'), l1 =3; end
    if strfind(e.roi2, 'ipsi'), l2 = 0;
    elseif strfind(e.roi2, 'contra'), l2 = 3; 
    elseif strfind(e.roi2, 'midline'), l2 = 3; end
    
    for ii=1:numel(freqs), if strcmp(freqs{ii}, e.freq), fidx = ii; end, end
    
    cohmean2{fidx}(l1+roi1, l2+roi2) = 100*mean(c.allcohC(:,idx)./c.allcohIC(:,idx)-1);
    cohstd2{fidx}(l1+roi1, l2+roi2) = 100*std(c.allcohC(:,idx)./c.allcohIC(:,idx)-1);
    ncomp2(l1+roi1, l2+roi2, fidx) = ncomp2(l1+roi1, l2+roi2, fidx)+1;
end
% average over connections within the same area
for k=1:numel(cohmean2)
    cohmean2{k} = cohmean2{k} + transpose(cohmean2{k});
    cohstd2{k} = cohstd2{k} + transpose(cohstd2{k});
    ncomp2(:,:,k) = ncomp2(:,:,k) + transpose(ncomp2(:,:,k));
end
    ncomp2(ncomp2(:)==0)=nan;
for k=1:numel(cohmean2)
   cohmean2{k} =  cohmean2{k}./ncomp2(:,:,k);
   cohmean2{k}(isnan(cohmean2{k})) = 0;
   cohstd2{k} =  cohstd2{k}./ncomp2(:,:,k);
   cohstd2{k}(isnan(cohstd2{k})) = 0;
end


maxnumber = abs(cat(3,cohmean2{:}));
maxnumber = max(maxnumber(:));
maxnumberstd = (cat(3,cohstd2{:}));
maxnumberstd = max(maxnumberstd(:));
for k=1:numel(cohmean2)
    figure;
    viscircles([0,0], 1, 'color','k'); hold on
    suptitle(sprintf('mean %s, nonsig',freqs{k}))
    vismot_circularGraph(cohmean2{k}(order, order),'std', cohstd2{k}(order, order), 'Colormap', repmat(cat(3,cmap(1,:), cmap(2,:)), [7 1 1]), 'Label', label, 'maxnumber', maxnumber, 'maxnumberstd', maxnumberstd);
end

%%%%%%%%%%%%
% T-VALUES %
%%%%%%%%%%%%
for k=1:numel(freqs), coh2{k} = zeros(7,7); end
ncomp=zeros(7,7,numel(freqs));
idx_nonsig = setdiff(1:numel(c.effect), c.idx_sign_uncorrected');
for k=1:numel(idx_nonsig)
    idx=idx_nonsig(k);
    e=c.effect(idx);
    
    % find ROI
    if strfind(e.roi1, 'occipital'), roi1 = 1;
    elseif strfind(e.roi1, 'parietal'), roi1 = 2;
    elseif strfind(e.roi1, 'motor'), roi1 = 3; 
    elseif strfind(e.roi1, 'frontal'), roi1 = 4; end
    if strfind(e.roi2, 'occipital'), roi2 = 1;
    elseif strfind(e.roi2, 'parietal'), roi2 = 2;
    elseif strfind(e.roi2, 'motor'), roi2 = 3; 
    elseif strfind(e.roi2, 'frontal'), roi2 = 4; end

    % find lateralization
    if strfind(e.roi1, 'ipsi'), l1 = 0;
    elseif strfind(e.roi1, 'contra'), l1 = 3; 
    elseif strfind(e.roi1, 'midline'), l1 =3; end
    if strfind(e.roi2, 'ipsi'), l2 = 0;
    elseif strfind(e.roi2, 'contra'), l2 = 3; 
    elseif strfind(e.roi2, 'midline'), l2 = 3; end
    
    for ii=1:numel(freqs), if strcmp(freqs{ii}, e.freq), fidx = ii; end, end
    
    coh2{fidx}(l1+roi1, l2+roi2) = e.stat + coh2{fidx}(l1+roi1, l2+roi2);
    ncomp(l1+roi1, l2+roi2, fidx) = ncomp(l1+roi1, l2+roi2, fidx)+1;
end
% average over connections within the same area
for k=1:numel(freqs)
    coh2{k} = coh2{k} + transpose(coh2{k});
    ncomp(:,:,k) = ncomp(:,:,k) + transpose(ncomp(:,:,k));
end
    ncomp(ncomp(:)==0)=nan;
for k=1:numel(freqs)
   coh2{k} =  coh2{k}./ncomp(:,:,k);
   coh2{k}(isnan(coh2{k})) = 0;
end

maxnumber = abs(cat(3,coh2{:}));
maxnumber = max(maxnumber(:));
for k=1:numel(freqs)
    figure;
    viscircles([0,0], 1, 'color','k'); hold on
    suptitle(sprintf('%s, not sig',freqs{k}))
    vismot_circularGraph(coh2{k}(order, order), 'Colormap', repmat(cat(3,cmap(1,:), cmap(2,:)), [7 1 1]), 'Label', label,'maxnumber', maxnumber);
end