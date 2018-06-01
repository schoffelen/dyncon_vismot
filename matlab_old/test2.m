pathname = '/raw/11/Project0012/JOE22/Seated1000/09-01-30@1414/';
runnames = {'2/';'3/';'4/';'5/';'6/'};

for k = 1:5
for m = 1:5
  fname = [pathname,runnames{k},'c,rfDC'];
  %event = read_event(fname);
  %hdr   = read_header(fname);
  
  cfg.datafile = fname;
  cfg.trialfun = ['trialfun_condition',num2str(m)];
  cd /analyse/4/Project0030/pilot
  [trlcfg{k,m}] = definetrial(cfg);
end
end

for m = 1:5
for k = 1:5
  cfg = trlcfg{k,m};
  cfg.channel = 'MEG';
  cfg.blc     = 'yes';
  if k==1,
    data = preprocessing(cfg);
  else
    data = appenddata([], data, preprocessing(cfg));
  end
end 

cfg              = [];
cfg.blc          = 'yes';
cfg.blcwindow    = [-0.5 0];
cfg.vartrllength = 2;
cfg.channel      = {'MEG' '-A40' '-A107' '-A157' '-A248'};
tlck{m}          = timelockanalysis(cfg, data);
end

%refchannel   = channelselection({'MEGREFL';'MEGREFG'}, hdr.label);
%cfg.channel  = refchannel;
%refdata      = preprocessing(cfg);

%cfg          = [];
%cfg.zscore   = 'yes';
%data         = denoise_pca(cfg,data,refdata);

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
lambda  = 0.001*s(1);

cfg = [];
cfg.method = 'lcmv';
cfg.grid   = grid;
cfg.vol    = vol;
cfg.fixedori = 'yes';
cfg.lambda = lambda;
cfg.keepfilter = 'yes';
source     = sourceanalysis(cfg, tlck2);
sd         = sourcedescriptives([], source);
