function [tlck,tlckall,tlckx] = doTimelockanalysis(subject);

cd(subject.pathname);
cd('data');
load([subject.name,'data']);

cfg              = [];
cfg.vartrllength = 2;
cfg.blc          = 'yes';
cfg.blcwindow    = [-0.2 0];
tlck(1)          = timelockanalysis(cfg, data1)
tlck(2)          = timelockanalysis(cfg, data2)
tlck(3)          = timelockanalysis(cfg, data3)
tlck(4)          = timelockanalysis(cfg, data4)
tlck(5)          = timelockanalysis(cfg, data5)

data    = appenddata([],data1,data2,data3,data4,data5);
data    = struct2double(data);
cfg.latency = [-0.2 0.6];
tlckall     = timelockanalysis(cfg,data);
cfg.keeptrials = 'yes';
tlckx       = timelockanalysis(cfg,data);
clear data;

%jackknife
strial = nansum(tlckx.trial,1);
ssmp   = reshape(tlckx.dof, [1 size(tlckx.dof)]);
for k = 1:size(tlckx.trial,1)
  tlckx.trial(k,:,:) = (strial-tlckx.trial(k,:,:))./(ssmp-1);
end


