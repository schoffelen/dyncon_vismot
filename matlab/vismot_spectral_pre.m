function [freqpre, tlckpre] = vismot_spectral_pre(subject,varargin)

% This function is based on doFreqanalysisMtmfft in the matlab_old folder
%
% do mtmfft analysis on data stratified for RT and balanced number of samples
% stimulus locked

smoothing = ft_getopt(varargin, 'smoothing', 4);
foilim    = ft_getopt(varargin, 'foilim', 'all');
output    = ft_getopt(varargin, 'output', 'pow');
doplanar  = ft_getopt(varargin, 'doplanar', 0);
conditions = ft_getopt(varargin, 'conditions', 'previous');
nrand      = ft_getopt(varargin, 'nrand', 100);

if smoothing==0
  taper = 'hanning';
  smoothing = 4;
else
  taper = 'dpss';
end

% load the data, this is the data ordered according to the conditions
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = vismot_data_reorder(alldata, conditions);

if ischar(foilim) && strcmp(foilim, 'all')
  foilim = [0 alldata.data1.fsample./2];
end

fd = fieldnames(alldata);
pre = cell(1,numel(fd));
for k = 1:numel(fd)
  data = alldata.(fd{k});
  
  cfg           = [];
  cfg.toilim    = [-0.5 0-1/256];
  cfg.minlength = 0.25;
  datapre       = ft_redefinetrial(cfg, data);
  clear data;
  
  cfg         = [];
  cfg.detrend = 'yes';
  datapre     = ft_preprocessing(cfg, datapre);
  pre{k}    = datapre;
  clear datapre;
end

if strcmp(output,'csd') 
	output = 'fourier';
	docsd  = true;
  dospectral = true;
elseif strcmp(output, 'tlck')
	docsd = false;
	dospectral = false;
else
	docsd = false;
	dospectral = true;
end

if dospectral
cfg         = [];
cfg.method  = 'mtmfft';
cfg.output  = output;
cfg.pad     = 600./pre{1}.fsample; %explicitly make nfft 600
cfg.foilim  = foilim;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
%cfg.channel = 'MEG';

% convert to synthetic planar gradient
if doplanar
  cfgn = [];
  cfgn.method   = 'template';
  cfgn.template = 'bti248_neighb.mat';
  neighbours = ft_prepare_neighbours(cfgn);
  
  cfgp = [];
  cfgp.method = 'sincos';
  cfgp.neighbours = neighbours;
  for k = 1:numel(pre)
    pre{k} = ft_megplanar(cfgp, pre{k});
  end
end

for k = 1:numel(pre)
  tmp = ft_freqanalysis(cfg, pre{k});
  if docsd
		ntap = sum(tmp.cumtapcnt);
		tmp  = ft_checkdata(tmp, 'cmbrepresentation', 'fullfast');
	  tmp.ntap = ntap;
	end
	freqpre(k) = tmp;
end

% combine planar gradients
if doplanar
  for k = 1:numel(pre)
    freqpre(k) = ft_combineplanar([], freqpre(k));
  end
end

else
  freqpre = [];	
end

% also compute the covariance matrix for lcmv
cfg = [];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
for k = 1:numel(pre)
	tlckpre(k) = ft_timelockanalysis(cfg, pre{k});
end


%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
