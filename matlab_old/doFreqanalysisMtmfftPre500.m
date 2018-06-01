function [freqpre] = doFreqanalysisiMtmfftPre500(subject,smoothing,foilim,flag);

%do mtmfft analysis on data stratified for RT and balanced number of samples
%stimulus locked

if smoothing==0,
  taper = 'hanning';
  smoothing = 4;
else
  taper = 'dpss';
end

if nargin<4,
  flag = 0;
end

cd(subject.pathname);
cd('data');
if flag==0,
  load([subject.name,'data']);
elseif flag==1,
  load([subject.name,'data_aligned']);
elseif flag==2,
  load([subject.name,'data_alignedclean']);
end

cd(subject.pathname);
cd('stratifyRT');
load([subject.name,'stratifyRT']);

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
  selrpt            = numel(input{k});
  trl               = findcfg(data.cfg, 'trl');
  trl(:,4)          = input{k}';
  data.cfg.trl      = trl;

  %cfg        = [];
  %cfg.trials = selrpt;
  %data       = preprocessing(cfg, data);

  if flag==2, data  = fixResample(data); end

  cfg           = [];
  cfg.toilim    = [-0.8 0-1/256];
  cfg.minlength = 0.8-2/256;
  data          = redefinetrial(cfg, data);
  cfg.toilim    = [-0.5 0-1/256];
  cfg.minlength = 0.5-2/256;
  data          = redefinetrial(cfg, data);

  cfg     = [];
  cfg.blc = 'yes';
  data    = preprocessing(cfg, data);

  Npre(k,1)     = length(data.trial);
  pre{k,1}      = data;
  clear data
end

cfg         = [];
cfg.method  = 'mtmfft';
cfg.output  = 'fourier';
%cfg.pad     = 256./pre{1}.fsample; %explicitly make nfft 256
%cfg.pad     = 'maxperlen'; %explicitly make nfft 256
cfg.pad     = 128./pre{1}.fsample;
cfg.foilim  = foilim;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
cfg.channel = 'MEG';
for k =1:4
  freqpre{k} = freqanalysis(cfg, pre{k});
end

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
