%cd /analyse/4/Project0030/freqstats;
cd /analyse/1/Project0002/tmpProject0030/;

d = dir;
d = d(3:end);
d = d(3:4:end); %only for tmpProject0030
names = {d(:).name}'
for k = 1:numel(names)
  names{k}
  load(names{k}, 'stat', 'stat13', 'stat42', 'statstrat');
  stat.powspctrm = stat.stat;
  stat13.powspctrm = stat13.stat;
  stat42.powspctrm = stat42.stat;
  statstrat.powspctrm = statstrat.stat;
  s{k,1} = stat;
  s13{k,1} = stat13;
  s42{k,1} = stat42;
  sstrat{k,1} = statstrat;
  clear stat stat13 stat42 statstrat;
end
freq  = selectdata(s{:}, 'param', 'powspctrm');
freq13  = selectdata(s13{:}, 'param', 'powspctrm');
freq42  = selectdata(s42{:}, 'param', 'powspctrm');
freqstrat = selectdata(sstrat{:}, 'param', 'powspctrm');
nsubj = size(freq.powspctrm,1);
freq.powspctrm(nsubj+[1:nsubj],:,:,:)   = 0;
freq13.powspctrm(nsubj+[1:nsubj],:,:,:) = 0;
freq42.powspctrm(nsubj+[1:nsubj],:,:,:) = 0;
freqstrat.powspctrm(nsubj+[1:nsubj],:,:,:) = 0;

design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];

cfg           = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'pooledTnan';
cfg.design    = design;
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.numrandomization = 0;
stat = ft_freqstatistics(cfg, freq);
stat13 = ft_freqstatistics(cfg, freq13);
stat42 = ft_freqstatistics(cfg, freq42);
statstrat = ft_freqstatistics(cfg, freqstrat);

