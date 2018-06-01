function [average, conditions] = doFreqstatsPlanarTFRperiRelchange(subject, bas)

%do freqbaseline of all conditions concatenated, and for each of the 
%conditions separately

fieldtripdefs

%cd([subject.pathname,'freq']);
%load([subject.name,'tfrperiHanning']);
cd('/analyse/1/Project0002/tmpProject0030');
load([subject.name,'tfrperiHanning2']);

cfg              = [];
cfg.baseline     = bas;
cfg.baselinetype = 'relchange';
for k = 1:numel(allfreq)
  tmp = selectdata(allfreq{k}, 'avgoverrpt', 'yes');
  tmp = ft_freqbaseline(cfg, tmp);
  conditions{k,1} = tmp;
  clear tmp;
end
tmp = selectdata(allfreq{:}, 'param', 'powspctrm', 'avgoverrpt', 'yes');
tmp = ft_freqbaseline(cfg, tmp);
average = tmp;
clear tmp;
