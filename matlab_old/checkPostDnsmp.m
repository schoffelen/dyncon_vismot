function [stat,freq,freq2] = checkPostDnsmp(subject,smoothing,foilim,conditions);

if smoothing==0,
  taper = 'hanning';
  smoothing = 4;
else
  taper = 'dpss';
end

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

cnt = 0;
for k = conditions
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
  clear data;

  cfg         = [];
  cfg.detrend = 'yes';
  datapst     = preprocessing(cfg, datapst);
  Npst = length(datapst.trial);

  cnt = cnt+1;
  pst{cnt} = datapst;
end

nsmp1 = cellfun('size',pst{1}.trial,2);
nsmp2 = cellfun('size',pst{2}.trial,2);

[ssmp1,ind1] = sort(nsmp1);
[ssmp2,ind2] = sort(nsmp2);

n = min(numel(ssmp1), numel(ssmp2));

[m,imx] = max([mean(ssmp1),mean(ssmp2)]);
if imx==1,
  indlong = ind1(end-n+1:end);
  smplong = ssmp1(end-n+1:end);
  indshort = ind2(1:n);
  smpshort = ssmp2(1:n);
  data     = pst{1};
elseif imx==2,
  indlong = ind2(end-n+1:end);
  smplong = ssmp2(end-n+1:end);
  indshort = ind1(1:n);
  smpshort = ssmp1(1:n);
  data     = pst{2};
end

cfg         = [];
cfg.method  = 'mtmfft';
cfg.output  = 'fourier';
cfg.pad     = 128./datapst.fsample; %explicitly make nfft 128
cfg.foilim  = foilim;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
cfg.channel = 'MEG';
cfg.trials  = indlong;
freq        = freqanalysis(cfg,data);

data2 = data;
for k = 1:length(indlong)
  tmpsmp = min(smplong(k),smpshort(k));
  data2.trial{indlong(k)} = blc(data.trial{indlong(k)}(:,1:tmpsmp));
  data2.time{indlong(k)} = data.time{indlong(k)}(1:tmpsmp);
end
freq2 = freqanalysis(cfg, data2);

freqx  = checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
freqx2 = checkdata(freq2, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.numrandomization = 0;
stat = freqstatistics(cfg, freqx, freqx2);
