%cd /analyse/4/Project0030/freqstats;
cd /analyse/1/Project0002/tmpProject0030/;

d = dir;
d = d(3:end);
d = d(1:2:end); %only for tmpProject0030
names = {d(:).name}'
for k = 1:numel(names)
  names{k}
  load(names{k}, 'avg');
  s{k,1} = avg;
  clear avg;
end
freq  = selectdata(s{:}, 'param', 'powspctrm');
nsubj = size(freq.powspctrm,1);
freq.powspctrm(nsubj+[1:nsubj],:,:,:) = 0;

design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];

cfg           = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.design    = design;
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.numrandomization = 0;
stat = ft_freqstatistics(cfg, freq);


