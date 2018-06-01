function [freq] = doFreqanalysisTFR(subject, frequency, tfwindow, smoothing, flag)

%compute high resolution tfr for a single frequency according to the additional
%input arguments

if nargin<5, flag = 0; end

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

if length(subject.runnames)>1,
  for k = 1:length(subject.runnames)
    fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{k},subject.datafile];
    event{k} = read_event(fname);
  end
else
  event{1} = read_event(fname);
end

cd(subject.pathname);
cd('data');
if flag,
  load([subject.name,'data_aligned']);
else
  load([subject.name,'data']);
end

if smoothing(1)==0,
  taper = 'hanning';
else
  taper = 'dpss';
end

for k = 1:5
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
  
  if k~=5,
    trl = findcfg(data.cfg, 'trl');
    if length(subject.runnames)>1,
      %make cell-array out of trl corresponding with the runs
      tmp  = diff([trl(:,1);0]);
      btrl = [0;find(tmp(1:end-1)<0)]+1;
      etrl = [find(tmp<0)];
      for kk = 1:numel(btrl)
        newtrl{kk,1} = trl(btrl(kk):etrl(kk),:);
      end
      trl = newtrl;
    else
      trl   = {trl};
    end

    rt = [];
    for kk = 1:length(event)
      cfg = [];
      cfg.eventtype   = 'RESPONSE';
      cfg.eventvalue  = 2.^(3*mod(k+1,2)+5);
      cfg.searchrange = 'afterzero';
      cfg.output      = 'samplefromoffset';
      tmprt = ft_recodeevent(cfg, event{kk}, trl{kk});
      rt    = [rt; tmprt(:)]; %in number of samples
    end
  else
    rt = nan+zeros(numel(data.trial),1); 
  end

  cfg           = [];
  cfg.method    = 'mtmconvol';
  cfg.output    = 'fourier';
  cfg.channel   = 'MEG';
  cfg.t_ftimwin = tfwindow;
  cfg.foi       = frequency;
  cfg.tapsmofrq = smoothing;
  cfg.taper     = taper;
  cfg.toi       = [-256:16:320]./data.fsample;
  cfg.pad       = (4*256)./data.fsample;
  tmpfreq       = ft_freqanalysis(cfg, data);
  tmp           = squeeze(sum(isnan(tmpfreq.fourierspctrm(:,1,1,:)),1));
  toiok         = find(tmp<size(tmpfreq.fourierspctrm,1));
  tmpfreq       = selectdata(tmpfreq, 'toilim', tmpfreq.time([toiok(1) toiok(end)]));
  tmpfreq.rt    = rt;
  freq{k,1}     = tmpfreq;
  clear tmpfreq;
end

function [x] = orthogonalise(input)

input = standardise(input,1);
x     = input(:,1);

for i = 2:size(input,2)
  y   = input(:,i);
  y   = y-x*(pinv(x'*x)*x'*y);
  if any(y), x = [x y]; end
end
