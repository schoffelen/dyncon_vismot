function [stat] = doFreqanalysisPlanarGLM(subject,smoothing,foilim);

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

  %match the pre and pst trials so that data can be stratified based on 
  %baseline power
  trlold  = datapre.cfg.trlold;
  trl     = datapre.cfg.trl;
  origvec = zeros(size(trlold,1),2);
  runold  = [0;find(diff(trlold(:,1))<0)]+1;
  runold(:,2) = [runold(2:end,1)-1;size(trlold,1)];
  run     = [0;find(diff(trl(:,1))<0)]+1;
  run(:,2) = [run(2:end,1)-1;size(trl,1)];
  for kk = 1:size(trl,1)
    block = find(kk<=run(:,2),1,'first');
    ii    = runold(block,1):runold(block,2);
    indx  = find(trlold(ii,1)<=trl(kk,1) & trlold(ii,2)>=trl(kk,2)) + runold(block,1) - 1;
    origvec(indx,1) = kk;
  end
  trlold  = datapst.cfg.trlold;
  trl     = datapst.cfg.trl;
  runold  = [0;find(diff(trlold(:,1))<0)]+1;
  runold(:,2) = [runold(2:end,1)-1;size(trlold,1)];
  run     = [0;find(diff(trl(:,1))<0)]+1;
  run(:,2) = [run(2:end,1)-1;size(trl,1)];
  for kk = 1:size(trl,1)
    block = find(kk<=run(:,2),1,'first');
    ii    = runold(block,1):runold(block,2);
    indx  = find(trlold(ii,1)<=trl(kk,1) & trlold(ii,2)>=trl(kk,2)) + runold(block,1) - 1;
    origvec(indx,2) = kk;
  end
  ok = origvec(:,1) > 0 & origvec(:,2) > 0;
  origvec = origvec(ok,:);

  cfg         = [];
  cfg.detrend = 'yes';
  cfg.trials  = origvec(:,1);
  datapre     = preprocessing(cfg, datapre);
  cfg.trials  = origvec(:,2);
  datapst     = preprocessing(cfg, datapst);
  Npre(k) = length(datapre.trial);
  Npst(k) = length(datapst.trial);

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

%npre = min(Npre);
%npst = min(Npst);
for k = 1:5
  %selvec = randperm(Npre(k));
  %selvec = selvec(1:npre);
  %Npre(k) = npre;
  %freqpre{k} = selectdata(freqpre{k}, 'rpt', selvec, 'channel', channel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', channel);
  %selvec = randperm(Npst(k));
  %selvec = selvec(1:npst);
  %Npst(k) = npst;
  %freqpst{k} = selectdata(freqpst{k}, 'rpt', selvec, 'channel', channel);
  freqpst{k} = selectdata(freqpst{k}, 'channel', channel);
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
cfg.method        = 'montecarlo';
cfg.statistic     = 'glm';
cfg.numrandomization = 0;
cfg.glm.statistic = 'T';
cfg.glm.contrast  = [1 0 0];

%condition 1: cue left, response left
%condition 2: cue left, response right
%condition 3: cue right, response left
%condition 4: cue right, response right

trlvec   = [0;cumsum([Npst(:);Npre(:)])];

tmpfd    = selectdata(fd, 'rpt', [trlvec(1)+1:trlvec(1+1) trlvec(3)+1:trlvec(3+1)]);
tmpfdbsl = selectdata(fd, 'rpt', [trlvec(6)+1:trlvec(6+1) trlvec(8)+1:trlvec(8+1)]);

design13(1,:) = [ones(1,Npst(1)) ones(1,Npst(3))*2];
design13(2,:) = freq.cumsumcnt([trlvec(1)+1:trlvec(1+1) trlvec(3)+1:trlvec(3+1)]);

nobs    = size(design13,2);
Tstat13 = zeros(length(tmpfd.label),length(tmpfd.freq));
for kk = 1:length(tmpfd.label)
  for mm = 1:length(tmpfd.freq)
    tmpdat = reshape(log10(tmpfd.powspctrm(:,kk,mm)),[1 nobs]);
    design13(3,:) = reshape(log10(tmpfdbsl.powspctrm(:,kk,mm)),[1 nobs]);
    tmp = statfun_glm(cfg, tmpdat, design13);
    Tstat13(kk,mm) = tmp.stat;
  end
end
 
tmpfd    = selectdata(fd, 'rpt', [trlvec(4)+1:trlvec(4+1) trlvec(2)+1:trlvec(2+1)]);
tmpfdbsl = selectdata(fd, 'rpt', [trlvec(9)+1:trlvec(9+1) trlvec(7)+1:trlvec(7+1)]);

design42(1,:) = [ones(1,Npst(4)) ones(1,Npst(2))*2];
design42(2,:) = freq.cumsumcnt([trlvec(4)+1:trlvec(4+1) trlvec(2)+1:trlvec(2+1)]);

nobs    = size(design42,2);
Tstat42 = zeros(length(tmpfd.label),length(tmpfd.freq));
for kk = 1:length(tmpfd.label)
  for mm = 1:length(tmpfd.freq)
    tmpdat = reshape(log10(tmpfd.powspctrm(:,kk,mm)),[1 nobs]);
    design42(3,:) = reshape(log10(tmpfdbsl.powspctrm(:,kk,mm)),[1 nobs]);
    tmp = statfun_glm(cfg, tmpdat, design42);
    Tstat42(kk,mm) = tmp.stat;
  end
end

stat = [];
stat.Tstat13 = Tstat13;
stat.Tstat42 = Tstat42;
stat.label   = tmpfd.label;
stat.freq    = tmpfd.freq;
stat.design13 = design13;
stat.design42 = design42;
stat.cfg      = [];
stat.dimord   = 'chan_freq';

