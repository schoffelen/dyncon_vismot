function [alltlck1, alltlck2] = doTimelockanalysisPlanar(subject, doplanar);

if nargin==1,
  doplanar = 1;
end

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

%get single trials for condition 1 and 3
cnt = 0;
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
  cnt               = cnt+1;
  data              = struct2double(data);

  cfg1 = [];
  cfg1.vartrllength = 2;
  cfg1.blc          = 'yes';
  cfg1.blcwindow    = [-0.5 0];
  cfg1.latency      = [-0.5 0.5-1/256];
  cfg1.channel      = 'MEG';
  cfg1.keeptrials   = 'yes';
  tlck1             = timelockanalysis(cfg1, data);

  if doplanar,
    cfg               = [];
    cfg.planarmethod  = 'sincos';
    cfg.neighbourdist = 0.037;
    data              = megplanar(cfg, data);
    tlck2             = timelockanalysis(cfg1, data);
  else
    tlck2 = [];
  end
  %cfg = [];
  %cfg.combinemethod = 'svd';
  %tlck              = combineplanar(cfg, tlck);

  alltlck1(cnt) = tlck1;
  if doplanar,
    alltlck2(cnt) = tlck2;
  else
    alltlck2(cnt) = nan;
  end
 
end

for k = 1:length(alltlck1)
  ind = sum(isfinite(alltlck1(k).trial(:,1,:)),3);
  sel = find(ind==256);
  alltlck1(k) = selectdata(alltlck1(k),'rpt',sel);
  if doplanar
    ind = sum(isfinite(alltlck2(k).trial(:,1,:)),3);
    sel = find(ind==256);
    alltlck2(k) = selectdata(alltlck2(k),'rpt',sel);
  end
end
