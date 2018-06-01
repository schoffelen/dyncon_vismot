function [stat,lrp,fd] = doFreqanalysisPlanar(subject,smoothing,foilim);

if smoothing==0,
  taper     = 'hanning';
  smoothing = 4;
else
  taper = 'dpss';
end

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

cnt = zeros(271,1);
for k = 1:5
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
  cfg.channel       = 'MEG';
  data              = preprocessing(cfg, data);
  label             = data.grad.label;
  [i1,i2]           = match_str(label,data.label);
  cnt(i1)           = cnt(i1) + 1;
  
  cfg               = [];
  cfg.planarmethod  = 'sincos';
  cfg.neighbourdist = 0.037;
  data              = megplanar(cfg, data);

  cfg           = [];
  cfg.toilim    = [0.2 0.7-1./data.fsample];
  cfg.minlength = 0.25;
  datapst       = redefinetrial(cfg, data);
  cfg.toilim    = [-0.5 0-1/data.fsample];
  datapre       = redefinetrial(cfg, data);
  clear data;

  cfg         = [];
  cfg.detrend = 'yes';
  datapre     = preprocessing(cfg, datapre);
  datapst     = preprocessing(cfg, datapst);
  Npre(k,1)   = length(datapre.trial);
  Npst(k,1)   = length(datapst.trial);

  cfg         = [];
  cfg.method  = 'mtmfft';
  cfg.output  = 'fourier';
  cfg.pad     = 128./datapre.fsample; %explicitly make nfft 128
  cfg.foilim  = foilim;
  cfg.taper   = taper;
  cfg.tapsmofrq = smoothing;
  freqpre{k}  = freqanalysis(cfg, datapre);
  freqpst{k}  = freqanalysis(cfg, datapst);
  clear datapre datapst
end
channel = label(cnt==5);
channel = [channel channel];
for k = 1:size(channel,1)
  channel{k,1} = [channel{k,1},'_dH'];
  channel{k,2} = [channel{k,2},'_dV'];
end
channel = channel(:);

npre = min(Npre);
npst = min(Npst);
for k = 1:5
  selvec = randperm(Npre(k));
  selvec = selvec(1:npre);
  Npre(k) = npre;
  freqpre{k} = selectdata(freqpre{k}, 'rpt', selvec, 'channel', channel);
  selvec = randperm(Npst(k));
  selvec = selvec(1:npst);
  Npst(k) = npst;
  freqpst{k} = selectdata(freqpst{k}, 'rpt', selvec, 'channel', channel);
end
freqpre = selectdata(freqpre{:}, 'param', 'fourierspctrm');
freqpst = selectdata(freqpst{:}, 'param', 'fourierspctrm');
freq    = selectdata(freqpst, freqpre, 'param', 'fourierspctrm');
clear freqpre freqpst;

cfg               = [];
cfg.combinemethod = 'svd';
freq              = combineplanar(cfg, freq);

cfg               = [];
cfg.channelcmb    = {};
cfg.keeptrials    = 'yes';
cfg.jackknife     = 'no';
fd                = freqdescriptives(cfg, freq);
 
cfg               = [];
cfg.design        = [ones(1,sum(Npst)) ones(1,sum(Npre))*2];
cfg.method        = 'montecarlo';
cfg.statistic     = 'indepsamplesT';
cfg.numrandomization = 0;
stat              = rmfield(freqstatistics(cfg, fd), 'cfg');

%condition 1: cue left, response left
%condition 2: cue left, response right
%condition 3: cue right, response left
%condition 4: cue right, response right

trlvec = [0;cumsum([Npst;Npre])];

tmpfd  = selectdata(fd, 'rpt', [trlvec(1)+1:trlvec(1+1) trlvec(3)+1:trlvec(3+1)]);
cfg.design  = [ones(1,Npst(1)) ones(1,Npst(3))*2];
stat(end+1) = rmfield(freqstatistics(cfg, tmpfd), 'cfg');
tmpfd  = selectdata(fd, 'rpt', [trlvec(4)+1:trlvec(4+1) trlvec(2)+1:trlvec(2+1)]);
cfg.design  = [ones(1,Npst(4)) ones(1,Npst(2))*2];
stat(end+1) = rmfield(freqstatistics(cfg, tmpfd), 'cfg');
%stat(end+1) = freqstatistics(cfg, tmpfd);

tmpfd  = selectdata(fd, 'rpt', [trlvec(6)+1:trlvec(6+1) trlvec(8)+1:trlvec(8+1)]);
cfg.design  = [ones(1,Npre(1)) ones(1,Npre(3))*2];
stat(end+1) = rmfield(freqstatistics(cfg, tmpfd), 'cfg');
%stat(end+1) = freqstatistics(cfg, tmpfd);
tmpfd  = selectdata(fd, 'rpt', [trlvec(9)+1:trlvec(9+1) trlvec(7)+1:trlvec(7+1)]);
cfg.design  = [ones(1,Npre(4)) ones(1,Npre(2))*2];
stat(end+1) = rmfield(freqstatistics(cfg, tmpfd), 'cfg');
%stat(end+1) = freqstatistics(cfg, tmpfd);


% load /home/jan/projects/visuomotor/matlab/lrplist
% cfg = [];
% cfg.channelcmb = channelcmb;
% avgL = stat(5);
% avgR = stat(6);
% avgL.time=avgL.freq;
% avgR.time=avgR.freq;
% avgL.avg=avgL.stat;
% avgR.avg=avgR.stat;
% lrp=lateralizedfield(cfg,avgL,avgR);
lrp = [];

%---
%stat(1) stimulus vs baseline all conditions collapsed
%stat(2) congruency effect response left
%stat(3) congruency effect response right
%stat(4) congruency effect response left baseline
%stat(5) congruency effect response right baseline
