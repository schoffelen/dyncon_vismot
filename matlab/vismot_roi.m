
foi = [12,22,46, 85];
datadir = '/project/3011085.03/analysis/source/mve/';
load atlas_subparc374_8k.mat
load cortex_inflated_shifted.mat;
atlas.pos = ctx.pos;


n=19;
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization = 2000;%5000;
cfgs.correctm = 'no';
cfgs.alpha = 0.025;
cfgs.design = [ones(1,n) 2*ones(1,n);1:n 1:n];
cfgs.correctm = 'no';
% cfgs.correctm='cluster';
% cfgs.clusterthreshold='nonparametric_individual';
% cfgs.connectivity = parcellation2connectivity_midline(atlas);
% cfgs.neighbours = cfgs.connectivity; %MvE
% cfgs.correcttail = 'prob';
% cfgs.clusteralpha = 0.05;
cfgs.parameter='stat';

cfgp=[];
cfgp.method = 'surface';
cfgp.funparameter = 'stat';
cfgp.funcolormap = flipud(brewermap(64, 'RdBu'));
cfgp.funcolorlim = 'maxabs';
cfgp.maskstyle = 'colormix';


cfg=[];
cfg.method = 'mean';
cfg.parameter = 'stat';

m=1;
for f=foi
    
suffix=sprintf('post_%d', f);
d   = dir(datadir);
sel = ~cellfun('isempty', strfind({d.name}', suffix));
d   = d(sel);

for k=1:numel(d)
    load(fullfile([datadir d(k).name]), 'statResp');
    resp{k} = statResp;
    resp{k}.pos  = atlas.pos;
    resp{k} = ft_sourceparcellate(cfg, resp{k}, atlas);
    nul{k} = resp{k};
    nul{k}.stat(:) = 0;
end

stat = ft_freqstatistics(cfgs, resp{:},nul{:});
stat.brainordinate = atlas;
source_parc{m} = stat;

ft_sourceplot(cfgp, source_parc{m}); material dull; h = light('position', [-1 0 -0.1]);
idx_max(m,:) = [0,0];
[val_max(m), idx_max(m,1)] = max(abs(source_parc{m}.stat(1:187,1))+abs(source_parc{m}.stat(188:end,1)));
idx_max(m,2) = idx_max(m,1)+187;
maxchan(m,:) = source_parc{m}.label(idx_max(m,:));
maxstat(m,:) = source_parc{m}.stat(idx_max(m,:));

showparcels = source_parc{1};
showparcels.stat(:)=0;
showparcels.stat(idx_max(m,:)) = 1;
% ft_sourceplot(cfgp, showparcels); material dull; h = light('position', [-1 0 -0.1]);

sprintf('%d Hz, max %s t=%d %d', f, maxchan{m}, maxstat(m,:))
% close all
m=m+1;
end
clear stat nul

s(1:4)=source_parc;

for k=1:4
    s{k}=s{k};
    s{k}.stat = (source_parc{k}.stat-nanmean(source_parc{k}.stat))./nanstd(source_parc{k}.stat);
end

Smot = s{1};
Smot.stat = (s{1}.stat+s{2}.stat-s{4}.stat)/3;
Svis = s{1};
Svis.stat = (-s{1}.stat+s{3}.stat+s{4}.stat)/3;

ft_sourceplot(cfgp, Smot); material dull; h = light('position', [-1 0 -0.1]);
ft_sourceplot(cfgp, Svis); material dull; h = light('position', [-1 0 -0.1]);

% combine left and right hemisphere to find symmetrical ROI's
Smot2=Smot;
Smot2.stat(1:187) = Smot.stat(1:187) - Smot.stat(188:end);
Smot2.stat(188:end) = Smot.stat(188:end) - Smot.stat(1:187);

Svis2=Svis;
Svis2.stat(1:187) = Svis.stat(1:187) - Svis.stat(188:end);
Svis2.stat(188:end) = Svis.stat(188:end) - Svis.stat(1:187);

% use ft_sourcemovie to find ROI's
cfgp=rmfield(cfgp, {'funcolorlim', 'maskstyle'});

Smot2 = rmfield(Smot2, 'freq');
Smot2.stat = [Smot2.stat Smot2.stat];
Smot2.time = 1:2;
Smot2.dimord = 'chan_time';
ft_sourcemovie(cfgp, Smot2);
idxmot = match_str(Smot.label, {'L_1_B05_02', 'L_1_B05_05', 'L_3_B05_02', 'L_4_B05_03', 'L_4_B05_07', 'L_2_B05_01', 'L_2_B05_04', 'R_1_B05_02', 'R_1_B05_05', 'R_3_B05_02', 'R_4_B05_03', 'R_4_B05_07', 'R_2_B05_01', 'R_2_B05_04'});
Smot.stat(setdiff(1:374,idxmot),1)=0;
ft_sourceplot(cfgp, Smot); material dull; h = light('position', [-1 0 -0.1]);


Svis2 = rmfield(Svis2, 'freq');
Svis2.stat = [Svis2.stat Svis2.stat];
Svis2.time = 1:2;
Svis2.dimord = 'chan_time';
ft_sourcemovie(cfgp, Svis2);
idxvis = match_str(Svis.label, {'R_19_B05_02','R_19_B05_04','R_19_B05_07','R_19_B05_08','R_19_B05_12','R_19_B05_14','R_18_B05_03','R_18_B05_04','L_19_B05_02','L_19_B05_04','L_19_B05_07','L_19_B05_08','L_19_B05_12','L_19_B05_14','L_18_B05_03','L_18_B05_04'});
Svis.stat(setdiff(1:374,idxvis),1)=0;
ft_sourceplot(cfgp, Svis); material dull; h = light('position', [-1 0 -0.1]);


