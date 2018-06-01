function [allfreq] = doFreqanalysisPlanar(subject,frequency,smoothing,twindow);

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

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
  cfg               = [];
  cfg.planarmethod  = 'sincos';
  cfg.neighbourdist = 0.037;
  data              = megplanar(cfg, data);

  cfg         = [];
  cfg.method  = 'mtmconvol';
  cfg.output = 'fourier';
  cfg.pad     = 1024./data.fsample; %explicitly make nfft 1024
  cfg.foi     = frequency;
  cfg.toi     = [fliplr([0:-4./data.fsample:-1]) 4./data.fsample:4./data.fsample:0.8];
  cfg.tapsmofrq = smoothing;
  cfg.t_ftimwin = twindow;
  cfg.channel   = 'MEG';
  freq          = freqanalysis(cfg, data);
  clear data

  cfg               = [];
  cfg.combinemethod = 'svd';
  freq              = combineplanar(cfg, freq);

  cfg               = [];
  cfg.channelcmb    = {};
  cfg.jackknife     = 'yes';
  fd                = freqdescriptives(cfg, freq);
  fd.nobs           = freq.cumtapcnt(1)*2*fd.dof-1; %the extra one is subtracted in the function
  fd                = fdsem2fdT(fd,[],'powspctrm',0,-0.25,'equalvar');

  allfreq(k)        = fd;
end

%condition 1: cue left, response left
%condition 2: cue left, response right
%condition 3: cue right, response left
%condition 4: cue right, response right
