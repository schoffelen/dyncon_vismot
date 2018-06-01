cd('/raw/11/Project0012/mss15/Pilot3/09-01-16@1639/2');
fname = [pwd,'/c,rfDC'];
event = read_event(fname);
hdr   = read_header(fname);

cfg.datafile = fname;
cfg.trialfun = 'mytrialfun';
cd /analyse/4/Project0030/pilot
cfg          = definetrial(cfg);
cfg.channel  = 'MEG';
cfg.blc      = 'yes';
data         = preprocessing(cfg);

refchannel   = channelselection({'MEGREFL';'MEGREFG'}, hdr.label);
cfg.channel  = refchannel;
refdata      = preprocessing(cfg);

cfg          = [];
cfg.zscore   = 'yes';
data         = denoise_pca(cfg,data,refdata);

%extract some channels which show blinks, and average
cfg          = [];
cfg.channel  = {'A194' 'A195' 'A228' 'A247'};
cfg.boxcar   = 0.1;
dataeog      = preprocessing(cfg, data);
dataeog      = checkdata(dataeog, 'datatype', 'timelock');
eog          = squeeze(mean(dataeog.trial,2));
eog          = reshape(standardise(eog(:),1),size(eog));
nsupra       = sum(eog>4,2);

cfg          = [];
cfg.trials   = find(nsupra==0);
cfg.dftfilter = 'yes';
data         = preprocessing(cfg, data);

cfg          = [];
cfg.method   = 'summary';
data         = rejectvisual(cfg, data);

cfg          = [];
cfg.detrend  = 'yes';
tlck         = timelockanalysis(cfg, data);
cfg            = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [0.20 0.25];
cfg.latency          = [0.20 0.25];
cfg.removemean       = 'no';
tlck2                = timelockanalysis(cfg, tlck);
cfg.covariancewindow = [-0.2 0];
cfg.latency          = [-0.2 0];
tlck3                = timelockanalysis(cfg, tlck);


load('mss15_grid6mm');
load('mss15_vol');

[a,b] = match_str(tlck2.label, tlck2.grad.label);
for k = 1:length(grid.inside)
  indx = grid.inside(k);
  grid.leadfield{indx} = grid.leadfield{indx}(b,:);
end

%compute a regularization parameter
[u,s,v] = svd(tlck2.cov);
lambda  = 0.01*s(1);

cfg = [];
cfg.method = 'lcmv';
cfg.grid   = grid;
cfg.vol    = vol;
cfg.fixedori = 'yes';
cfg.lambda = '10%';
cfg.keepfilter = 'yes';
source     = sourceanalysis(cfg, tlck2);
sd         = sourcedescriptives([], source);
