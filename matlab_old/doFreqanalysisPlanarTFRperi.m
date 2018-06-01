function [allfreq] = doFreqanalysisPlanarTFRperi(subject,frequency,smoothing,twindow,toi);

if nargin==4,
  toi     = [-0.79:0.01:0.78];
end

cd(subject.pathname);
cd('data');
load([subject.name,'data']);
cd(subject.pathname);
cd('rt');
load([subject.name,'rt']);

if all(smoothing==0),
  taper = 'hanning';
else
  taper = 'dpss';
end  

for k = 1:4
  warning off;
  if k==1,
    data = data1;
  elseif k==2,
    data = data2;
  elseif k==3,
    data = data3;
  elseif k==4,
    data = data4;
  elseif k==5,
    data = data5;
  end
  data              = struct2double(data);
  %remove rounding error from time axes introduced by resampling
  for m = 1:numel(data.trial)
    data.time{m} = round(256*data.time{m})./256;
  end
  trl = findcfg(data.cfg, 'trl');
  trl(:,4) = rt{k};
  data.cfg.trl = trl;

  cfg        = [];
  cfg.toilim = [-1 1-1./256];
  data       = ft_redefinetrial(cfg, data);

  cfg               = [];
  cfg.planarmethod  = 'sincos';
  cfg.neighbourdist = 0.037;
  data              = ft_megplanar(cfg, data);

  cfg         = [];
  cfg.detrend = 'yes';
  data        = ft_preprocessing(cfg, data);

  cfg         = [];
  cfg.method  = 'mtmconvol';
  cfg.output = 'pow';
  cfg.pad     = 2; %explicitly make nfft 1024
  cfg.foi     = frequency;
  cfg.toi     = toi;
  cfg.tapsmofrq = smoothing;
  cfg.t_ftimwin = twindow;
  cfg.taper     = taper;
  cfg.precision = 'single';
  cfg.channel   = 'MEG';
  cfg.keeptrials = 'yes';
  freq          = ft_freqanalysis(cfg, data);
  clear data

  cfg               = [];
  freq              = combineplanar(cfg, freq);

  allfreq{k}        = freq;
end

%condition 1: cue left, response left
%condition 2: cue left, response right
%condition 3: cue right, response left
%condition 4: cue right, response right
