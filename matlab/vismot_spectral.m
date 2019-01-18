function [freq, tlck] = vismot_spectral(subject,varargin)

% This function is based on doFreqanalysisMtmfft in the matlab_old folder
%
% do mtmfft analysis on data stratified for RT and balanced number of samples
% stimulus locked
toi        = ft_getopt(varargin, 'toi', 'post');
smoothing  = ft_getopt(varargin, 'smoothing', 4);
foilim     = ft_getopt(varargin, 'foilim', 'all');
output     = ft_getopt(varargin, 'output', 'pow');
doplanar   = ft_getopt(varargin, 'doplanar', 0);
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
data_short = cell(1,numel(fd));
if strcmp(toi, 'post')
    toilim = [0 0.5-1/256];
elseif strcmp(toi, 'pre')
    toilim = [-0.5 0-1/256];
else
    error('please specificy toi as *pre* or *post')
end
for k = 1:numel(fd)
  data = alldata.(fd{k});
  
  cfg           = [];
  cfg.toilim    = toilim;
  cfg.minlength = 0.5;
  data_tmp    = ft_redefinetrial(cfg, data);
  clear data;
  
  cfg         = [];
  cfg.detrend = 'yes';
  data_tmp     = ft_preprocessing(cfg, data_tmp);
  data_short{k}    = data_tmp;
  clear data_tmp;
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
cfg.pad     = 600./data_short{1}.fsample; %explicitly make nfft 600
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
  for k = 1:numel(data_short)
    data_short{k} = ft_megplanar(cfgp, data_short{k});
  end
end

for k = 1:numel(data_short)
  tmp = ft_freqanalysis(cfg, data_short{k});
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
