function [stat13, stat42, stat, stat13strat, stat42strat, statstrat] = doFreqstatsPlanarTFRperiCong(subject, bas)

%do freqbaseline of all conditions concatenated, and for each of the 
%conditions separately

fieldtripdefs

%cd([subject.pathname,'freq']);
%load([subject.name,'tfrperiHanning']);
cd('/analyse/1/Project0002/tmpProject0030/');
load([subject.name,'tfrperiHanning2']);

cfg              = [];
cfg.baseline     = bas;
cfg.baselinetype = 'relchange';
for k = 1:numel(allfreq)
  allfreq{k} = ft_freqbaseline(cfg, allfreq{k});
end

cfg        = [];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.parameter = 'powspctrm';
cfg.numrandomization = 0;
stat42 = ft_freqstatistics(cfg,allfreq{4},allfreq{2});
stat13 = ft_freqstatistics(cfg,allfreq{1},allfreq{3});

load('/home/jan/projects/visuomotor/matlab/lrplist');
[a1,b1] = match_str(channelcmb(:,1), stat42.label);
[a2,b2] = match_str(channelcmb(:,2), stat13.label);
stat    = stat42;
stat.label = channelcmb(:,1);
stat.stat  = zeros(numel(stat.label),numel(stat.freq),numel(stat.time));
stat.stat  = (stat42.stat(b1,:,:) + stat13.stat(b2,:,:))./sqrt(2);

[input, output, binaxis] = stratifyRT(subject, allfreq);
for k = 1:4
  allfreq{k} = selectdata(allfreq{k},'rpt',find(isfinite(output{k})));
end
stat42strat = ft_freqstatistics(cfg,allfreq{4},allfreq{2});
stat13strat = ft_freqstatistics(cfg,allfreq{1},allfreq{3});
statstrat   = stat42strat;
statstrat.label = channelcmb(:,1);
statstrat.stat  = zeros(numel(statstrat.label),numel(statstrat.freq),numel(statstrat.time));
statstrat.stat  = (stat42strat.stat(b1,:,:) + stat13strat.stat(b2,:,:))./sqrt(2);
