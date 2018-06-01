function [freq, freqbeam] = doTFR(subject, frequency, covariancewindow)

cd([subject.pathname,'data']);
try,
  load([subject.name,'data_aligned']);
catch
  load([subject.name]);
end

for k = 1:4
  warning off;
  if k==1,
    data = struct2double(data1);
  elseif k==2,
    data = struct2double(data2);
  elseif k==3,
    data = struct2double(data3);
  elseif k==4,
    data = struct2double(data4);
  elseif k==5,
    data = struct2double(data5);
  end
  warning on;
  
  %round the time points
  for m = 1:numel(data.trial)
    data.time{m} = round(data.time{m}*data.fsample)./data.fsample;
  end

  cfg     = [];
  cfg.method = 'mtmconvol';
  cfg.output = 'fourier';
  cfg.channel = ft_channelselection({'MEG'}, data.label);
  cfg.toi = [-128:16:256]./data.fsample;
  cfg.foi = frequency(1);
  cfg.t_ftimwin = frequency(2);
  cfg.tapsmofrq = frequency(3);
  if cfg.tapsmofrq==0
    cfg.taper     = 'hanning';
    cfg.tapsmofrq = 1./cfg.t_ftimwin;
  end
  cfg.pad =  4*256./data.fsample;
  freq = ft_freqanalysis(cfg, data);
end
