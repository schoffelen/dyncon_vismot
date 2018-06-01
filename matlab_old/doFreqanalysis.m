function [freqpre,freqpst] = doFreqanalysis(subject,smoothing,foilim,flag);

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

  cfg           = [];
  cfg.toilim    = [0.2 0.7-1./data.fsample];
  cfg.minlength = 0.25;
  datapst       = redefinetrial(cfg, data);
  cfg.toilim    = [-0.5 0-1/256];
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
  Npre = length(datapre.trial);
  Npst = length(datapst.trial);

  cfg         = [];
  cfg.method  = 'mtmfft';
  cfg.output  = 'fourier';
  cfg.pad     = 128./datapre.fsample; %explicitly make nfft 128
  cfg.foilim  = foilim;
  cfg.taper   = taper;
  cfg.tapsmofrq = smoothing;
  cfg.channel = 'MEG';
  pre     = freqanalysis(cfg, datapre);
  pst     = freqanalysis(cfg, datapst);
  clear datapre datapst
  
  %pre = selectdata(pre, 'foilim', [10 20 56]);
  %pst = selectdata(pst, 'foilim', [10 20 56]);

  %cfg            = [];
  %cfg.channelcmb = {};
  %fdpre          = freqdescriptives(cfg, freqpre);
  %fdpst          = freqdescriptives(cfg, freqpst);
  %fdpre.nobs     = 2*sum(freqpre.cumtapcnt)-2;
  %fdpst.nobs     = 2*sum(freqpst.cumtapcnt)-2;
  %fdx            = fdsem2fdT(fdpst,fdpre,'powspctrm',0,[],'equalvar');

  freqpst{k} = pst;
  freqpre{k} = pre;

end

%condition 1: cue left, response left
%condition 2: cue left, response right
%condition 3: cue right, response left
%condition 4: cue right, response right
