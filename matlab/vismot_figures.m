addpath('/project/3011085.03/scripts/RainCloudPlots/tutorial_matlab/')
addpath('/project/3011085.03/scripts/Robust_Statistical_Toolbox/')
subjects = vismot_subjinfo;

filename = [subjects(1).pathname, '/rt/', 'stat_simon.mat'];
d = load(filename);

cmap = flipud(brewermap(2,'RdBu'));
f1 = figure;
h1 = raincloud_plot(d.avgC, 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);
h2 = raincloud_plot(d.avgIC, 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
legend([h1{1} h2{1}], {'Group 1', 'Group 2'});
title(['Figure 7' newline 'A) Dodge Options Example 1']);
set(gca,'XLim', [0 40], 'YLim', [-.075 .15]);
box off

cmap = brewermap(2, 'RdBu');
filename = [subjects(1).pathname, '/rt/', 'stat_gratton.mat'];
d2 = load(filename);
figure;
data{1,1} = d2.avgC_C;
data{2,1} = d2.avgN_C;
data{3,1} = d2.avgIC_C;
data{1,2} = d2.avgC_IC;
data{2,2} = d2.avgN_IC;
data{3,2} = d2.avgIC_IC;
h  = rm_raincloud(data, cmap);
% note that the x and y axes are switched
xlabel('Reaction time (s)')
yticklabels({'previous IC', 'previous N', 'previous C'}) % in reverse order because of switched axes!
xlim([400 1100])
legend({'current congruent', 'current incongruent'}, 'Location', 'NorthEast') % how to implement this?



%% Post cue power
alldir = '/project/3011085.03/';
load(fullfile([alldir, 'analysis/source/roi.mat']));
load(fullfile([alldir, 'analysis/stat_bf_post.mat']))

% make slice plots of power


% Make bar graph of defined of power in defined ROIs and FOIs
% find indices of roi locations
for k=1:size(roi,1)-1
    idx_left(k) = find_dipoleindex(stat, roi{k+1,3});
    idx_right(k) = find_dipoleindex(stat, roi{k+1,4});
end

cmap = brewermap(5, 'Set3');
figure;
% occipital
subplot(1,3,1);
occ = [stat.stat(idx_left(1),1), stat.stat(idx_left(2),3) stat.stat(idx_left(3),4), stat.stat(idx_left(4),5); stat.stat(idx_right(1),1), stat.stat(idx_right(2),3) stat.stat(idx_right(3),4), stat.stat(idx_right(4),5)];
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
par = [stat.stat(idx_left(5),1), stat.stat(idx_left(6),3) stat.stat(idx_left(7),4), stat.stat(idx_left(8),5); stat.stat(idx_right(5),1), stat.stat(idx_right(6),3) stat.stat(idx_right(7),4), stat.stat(idx_right(8),5)];
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
mot = [stat.stat(idx_left(9),1), stat.stat(idx_left(10),2) stat.stat(idx_left(11),3), stat.stat(idx_left(12),5); stat.stat(idx_right(9),1), stat.stat(idx_right(10),2) stat.stat(idx_right(11),3), stat.stat(idx_right(12),5)];
b = bar(mot);
tmpcmap = cmap([1:3 5],:);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'beta', 'gamma1', 'gamma3'}, 'location', 'northeast')
title('motor')

%% ROI plot


%% pre cue power

s = load(fullfile([alldir, 'analysis/stat_bf_pre.mat']))
% OPTION 1
cmap = brewermap(12, 'Set1');
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
figure; iosr.statistics.boxPlot(diff', 'showScatter', true, 'scatterMarker', '.', 'boxColor'); 
load(fullfile([alldir, 'analysis/source/roi.mat']));
cmap = brewermap(5, 'Set3');
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
b = bar([lpow(1:4); rpow(1:4)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
subplot(1,3,2)
b = bar([lpow(5:8); rpow(5:8)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});
subplot(1,3,3)
tmpcmap = cmap([1:3 5],:);
b = bar([lpow(9:12); rpow(9:12)]);
for k=1:4
  set(b(k),'FaceColor', tmpcmap(k,:));
end
set(gca,'xticklabel',{'left', 'right'});