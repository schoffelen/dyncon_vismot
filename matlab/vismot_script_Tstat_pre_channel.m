datadir = '/project/3011085.03/analysis/freq/mve/';
foilim=[0 60];
% foilim=[40 120];
whichstat = 'statResp';
d = dir(fullfile(datadir,sprintf('*axial*%d-%d.mat',foilim(1), foilim(2))));

for m = 1:numel(d)
        dum = load(fullfile(d(m).folder,d(m).name),whichstat);
        tmp(m) = dum.(whichstat);
end
dat = permute(cat(3,tmp.stat), [3,1,2]);
n = size(dat,1);
for m = 1:n
        dat(m+n,:,:) = repmat(nanmean(nanmean(dat(m, :, :))), [size(dat,2), size(dat,3)]);
end

data.label = dum.(whichstat).label;
data.freq = dum.(whichstat).freq;
data.powspctrm = dat;
data.dimord = 'rpt_chan_freq';

cfg=[];
cfg.design   = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.numrandomization = 1000;
stat      = ft_freqstatistics(cfg, data);
    
cfgp=[];
cfgp.layout = '4D248_helmet.mat';
cmap = flipud(brewermap(64,'RdBu'));
cfgp.colormap = cmap;
cfgp.parameter = 'stat';
ft_singleplotER(cfgp, stat);
