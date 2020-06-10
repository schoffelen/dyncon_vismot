function [freq] = vismot_spectral_prepost_long(subject,varargin)

% do mtmconvol analysis on data 

range     = ft_getopt(varargin, 'range', 'low');
doplanar  = istrue(ft_getopt(varargin, 'doplanar', false));


% load the data
load(fullfile(subject.pathname,'data',[subject.name,'data_long']));

pre = cell(1,5);
pst = cell(1,5);
for k = 1:5
  if k==1
    data = data1;
  elseif k==2
    data = data2;
  elseif k==3
    data = data3;
  elseif k==4
    data = data4;
  elseif k==5
    data = data5;
  end
 
  cfg         = [];
  cfg.detrend = 'yes';
  data        = ft_preprocessing(cfg, data);
  dat{k}      = data;
  clear data;
end

switch range
  case 'low'
    foi = 2.5:2.5:40;
    taper = 'hanning';
    smoothing = ones(1,numel(foi));
    t_ftimwin = ones(1,numel(foi))./2.5;
  case 'high'
    foi = 30:4:100;
    taper = 'dpss';
    smoothing = ones(1,numel(foi)).*8;
    t_ftimwin = ones(1,numel(foi))./4;
end

cfg         = [];
cfg.method  = 'mtmconvol';
cfg.pad     = 1200./dat{1}.fsample; %explicitly make nfft 300
cfg.foi     = foi;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
cfg.channel   = dat{1}.label;
cfg.t_ftimwin = t_ftimwin;
cfg.keeptrials = 'yes';
cfg.toi    = (-150:15:180)./300;

% convert to synthetic planar gradient
if doplanar
  cfgn = [];
  cfgn.method   = 'template';
  cfgn.template = 'bti248_neighb.mat';
  neighbours = ft_prepare_neighbours(cfgn);
  
  cfgp = [];
  cfgp.method = 'sincos';
  cfgp.neighbours = neighbours;
  for k = 1:5
    dat{k} = ft_megplanar(cfgp, dat{k});
  end
end

for k = 1:5
  freq(k) = ft_freqanalysis(rmfield(cfg, 'channel'), dat{k});
end

% combine planar gradients
if doplanar
  for k = 1:5
    freq(k) = ft_combineplanar([], freq(k));
  end
end

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
