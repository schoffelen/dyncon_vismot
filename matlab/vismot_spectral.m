function [freq, tlck] = vismot_spectral(subject,varargin)

% This function is based on doFreqanalysisMtmfft in the matlab_old folder
%
% do mtmfft analysis on data stratified for RT and balanced number of samples
% stimulus locked: Note that this functionality had disappeared in the
% reorganization, so it hasn't been used for a while. It is back now (Feb
% 05, 2019)

toi        = ft_getopt(varargin, 'toi', 'post');
conditions = ft_getopt(varargin, 'conditions', []);
smoothing  = ft_getopt(varargin, 'smoothing', 4);
foilim     = ft_getopt(varargin, 'foilim', 'all');
output     = ft_getopt(varargin, 'output', 'pow');
doplanar    = istrue(ft_getopt(varargin, 'doplanar',  false));
doprewhiten = istrue(ft_getopt(varargin, 'prewhiten', false));

% this determines whether the trials are going to be divided into
% structures according to the condition of the current, or previous trial
if isempty(conditions)
  if strcmp(toi, 'post')
    conditions = 'current';
  elseif strcmp(toi, 'pre')
    conditions = 'previous';
  end
end

if smoothing==0
  taper = 'hanning';
else
  taper = 'dpss';
end

% load the data, this is the data ordered according to the condition of the
% current trial
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));

% reorder if needed
alldata = vismot_data_reorder(alldata, conditions);

if ischar(foilim) && strcmp(foilim, 'all')
  foilim = [0 alldata.data1.fsample./2];
end

if strcmp(toi, 'post')
  toilim = [0.2 0.7-1/300];
elseif strcmp(toi, 'pre')
  toilim = [-0.5 0-1/300]; % don't see a reason to use 1/256 if the resample fs is 300
elseif strcmp(toi, 'prepost')
  toilim = [0.2 0.7-1/300; -0.5 0-1/300];
else
  error('please specificy toi as *pre* or *post*, or *prepost*')
end

fd         = fieldnames(alldata);
data_short = cell(size(toilim,1),numel(fd));
for k = 1:numel(fd)
  data = alldata.(fd{k});
  
  for m = 1:size(toilim,1)
    cfg           = [];
    cfg.toilim    = toilim(m,:);
    cfg.minlength = 0.5;
    data          = ft_redefinetrial(cfg, data); % note: this should actually use ft_selectdata, but for some reason this does not work robustly, due to rounding issues of time axes or so
    
    cfg           = [];
    cfg.detrend   = 'yes';
    data          = ft_preprocessing(cfg, data);
    data_short{m,k} = data;
    clear data;
  end
end

if strcmp(output,'csd') 
	output = 'fourier';
	docsd  = true;
  dospectral = true;
elseif strcmp(output, 'tlck')
	docsd      = false;
	dospectral = false;
else
	docsd      = false;
	dospectral = true;
end

if dospectral
  cfg         = [];
  cfg.method  = 'mtmfft';
  cfg.output  = output;
  cfg.pad     = 600./data_short{1}.fsample; %explicitly make nfft 600
  cfg.foilim  = foilim;
  cfg.taper   = taper;
  cfg.tapsmofrq = smoothing;
  
  if doprewhiten
    load(fullfile(subject.pathname,'data',[subject.name,'emptyroom']));
    
    cfgr        = [];
    cfgr.length = 0.5;
    noise       = ft_redefinetrial(cfgr, data);
    clear data;
    
    cfgd         = [];
    cfgd.detrend = 'yes';
    noise        = ft_preprocessing(cfgd, noise);
  end
  
  if doplanar
    % convert to synthetic planar gradient
    cfgn          = [];
    cfgn.method   = 'template';
    cfgn.template = 'bti248_neighb.mat';
    neighbours    = ft_prepare_neighbours(cfgn);
    
    cfgp            = [];
    cfgp.method     = 'sincos';
    cfgp.neighbours = neighbours;
    for k = 1:numel(data_short)
      data_short{k} = ft_megplanar(cfgp, data_short{k});
    end
    if doprewhiten
      noise = ft_megplanar(cfgp, noise);
    end
  end
  
  if doprewhiten
    noise = ft_freqanalysis(cfg, noise);
    noise = ft_checkdata(noise, 'cmbrepresentation', 'fullfast');
  end

  for k = 1:numel(data_short)
    tmp = ft_freqanalysis(cfg, data_short{k});
    if doprewhiten
      tmp = ft_denoise_prewhiten([], tmp, noise);
    end
    if docsd
      ntap = sum(tmp.cumtapcnt);
      tmp  = ft_checkdata(tmp, 'cmbrepresentation', 'fullfast');
      tmp.ntap = ntap;
    end
    freq(k) = tmp;
  end

  % combine planar gradients
  if doplanar
    for k = 1:numel(data_short)
      freq(k) = ft_combineplanar([], freq(k));
    end
  end

else
  freq = [];	
end

% also compute the covariance matrix for lcmv
cfg = [];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
for k = 1:numel(data_short)
	tlck(k) = ft_timelockanalysis(cfg, data_short{k});
end


%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
