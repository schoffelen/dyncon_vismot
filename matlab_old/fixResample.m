function [data] = fixResample(data)

%function to convert old-style resampled data into new-style, i.e.
%constructing a cfg.resampletrl matrix and append it to the data's
%configuration

trl   = findcfg(data.cfg, 'trl');
if size(trl,1) ~= length(data.trial)
  trl = findcfg(data.cfg.previous{2}, 'trl');
end
fsold = findcfg(data.cfg, 'origfs');
fsnew = findcfg(data.cfg, 'resamplefs');
nsmp  = trl(:,2)-trl(:,1)+1;

for k = 1:size(trl,1)
  time{k} = offset2time(trl(k,3),fsold,nsmp(k));
end

dumdata = [];
dumdata.fsample = fsold;
dumdata.label   = 'time';
dumdata.time    = time;
dumdata.trial   = time;
dumdata.cfg.trl = trl;

cfg            = [];
cfg.resamplefs = fsnew;
cfg.detrend    = 'no';
cfg.blc        = 'yes';
dumdatax       = resampledata(cfg, dumdata);

data.cfg.resampletrl = findcfg(dumdatax.cfg, 'resampletrl');
