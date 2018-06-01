function [alltlck] = doTimelockanalysisPlanar(subject);

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

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
  cfg.planarmethod  = 'sincos';
  cfg.neighbourdist = 0.037;
  data              = megplanar(cfg, data);

  cfg = [];
  cfg.vartrllength = 2;
  cfg.blc          = 'yes';
  cfg.blcwindow    = [-0.5 0];
  cfg.latency      = [-0.5 0.5-1/256];
  cfg.channel       = 'MEG';
  tlck             = timelockanalysis(cfg, data);

  cfg = [];
  %cfg.combinemethod = 'svd';
  tlck              = combineplanar(cfg, tlck);

  alltlck(k) = tlck;
end

%allavg = cat(2,alltlck(:).avg);
%tlckx  = alltlck(1);
%tlckx.time = cat(2,alltlck(:).time);
%tlckx.avg  = allavg;
%cfg = [];
%cfg.combinemethod = 'svd';
%tlckx             = combineplanar(cfg, tlckx);
%
%alltlck(1).avg = tlckx.avg(:,1:length(alltlck(1).time));
%alltlck(1).label = tlckx.label;
%tlckx.avg(:,1:length(alltlck(1).time)) = [];
%alltlck(2).avg = tlckx.avg(:,1:length(alltlck(2).time));
%alltlck(2).label = tlckx.label;
%tlckx.avg(:,1:length(alltlck(2).time)) = [];
%alltlck(3).avg = tlckx.avg(:,1:length(alltlck(3).time));
%alltlck(3).label = tlckx.label;
%tlckx.avg(:,1:length(alltlck(3).time)) = [];
%alltlck(4).avg = tlckx.avg(:,1:length(alltlck(4).time));
%alltlck(4).label = tlckx.label;
%tlckx.avg(:,1:length(alltlck(4).time)) = [];
%
%
%alltlck(1).label = tlckx.label;
