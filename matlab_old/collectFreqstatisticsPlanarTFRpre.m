function [stat1, stat2, stat3] = collectFreqstatisticsPlanarTFRpre(pathname, suffix, nrand)

cd(pathname);
d = dir;
d = d(3:end);

names = {d(:).name}';
sel   = find(~cellfun('isempty', strfind(names, suffix)));
names = names(sel);

for k = 1:numel(names)
  names{k}
  load(names{k});
  stat13.powspctrm = stat13.stat; stat13 = rmfield(stat13, 'stat');
  stat42.powspctrm = stat42.stat; stat42 = rmfield(stat42, 'stat');
  s13{k,1} = stat13;
  s42{k,1} = stat42;
  clear stat13 stat42
end

stat13 = selectdata(s13{:}, 'param', 'powspctrm');
stat42 = selectdata(s42{:}, 'param', 'powspctrm');

%pool the response conditions
load lrplist
[a,b] = match_str(channelcmb(:,1), stat42.label);
stat  = stat42;
stat.powspctrm = stat.powspctrm(:,b,:,:);
stat.label     = stat.label(b);
[a,b] = match_str(channelcmb(:,2), stat13.label);
stat.powspctrm = (stat.powspctrm + stat13.powspctrm(:,b,:,:))./sqrt(2);

nsubj  = size(stat.powspctrm,1);
design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];
stat.powspctrm(nsubj+1:2*nsubj,:,:,:) = 0;

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.method    = 'montecarlo';
cfg.numrandomization = 1000;
cfg.statistic = 'pooledT';
cfg.correctm  = 'no';
cfg.design    = design;
cfg.ivar      = 1;
cfg.uvar      = 2;
sx            = freqstatistics(cfg, stat);
