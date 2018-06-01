function [stat] = doFreqanalysisPlanar(subject,smoothing);

cd(subject.pathname);
cd('freq');
load([subject.name,'mtmfft',num2str(smoothing,'%03d')]);
load([subject.name,'stratify',num2str(smoothing,'%03d')]);

%subselect channels
for k = 1:length(freqpre)
  freqpre{k} = selectdata(freqpre{k}, 'channel', label);
  freqpst{k} = selectdata(freqpst{k}, 'channel', label);
end

%convert to planar 
cfg               = [];
cfg.neighbourdist = 0.037;
cfg.planarmethod  = 'sincos';
for k = 1:length(freqpre)
  freqpre{k} = megplanar(cfg, freqpre{k});
  freqpst{k} = megplanar(cfg, freqpst{k});
end

freq = cat(2,freqpst, freqpre);
clear freqpst freqpre

for k = 1:length(freq)
  n(k,1) = length(freq{k}.cumtapcnt);
end
freq   = selectdata(freq{:}, 'param', 'fourierspctrm');

%%combine planar with a particular contrast in mind and do statistics
%cfg1               = [];
%cfg1.combinemethod = 'svd';
%freq               = combineplanar(cfg1, freq);

cfg2               = [];
cfg2.channelcmb    = {};
cfg2.keeptrials    = 'yes';
cfg2.jackknife     = 'no';
freq               = freqdescriptives(cfg2, freq);

cfg1               = [];
cfg1.combinemethod = 'sum';
freq               = combineplanar(cfg1, freq);

cfg3                  = [];
cfg3.method           = 'montecarlo';
cfg3.statistic        = 'indepsamplesT';
cfg3.numrandomization = 0;

alldesign = zeros(1,0);
for k = 1:length(n)
  alldesign(end+1:end+n(k)) = k;
end
contrasts = [1 3;4 2;6 8;9 7;1 6;2 7;3 8;4 9;5 10;1 5;2 5;3 5;4 5];
%1 = congruency effect left
%2 = congruency effect right
%3 = congruency effect baseline left
%4 = congruency effect baseline left
%5 = stimulus vs baseline condition 1
%6 = stimulus vs baseline condition 2
%7 = stimulus vs baseline condition 3
%8 = stimulus vs baseline condition 4
%9 = stimulus vs baseline condition 5
%10 = condition 1 vs condition 5
%11 = condition 2 vs condition 5
%12 = condition 3 vs condition 5
%13 = condition 4 vs condition 5
freq.powspctrm = log10(freq.powspctrm);
for cc = 1:size(contrasts,1)
  warning off
  
  %tmpfreq = freq(contrasts(cc,:));
  %n1      = length(tmpfreq{1}.cumtapcnt);
  %n2      = length(tmpfreq{2}.cumtapcnt);
  %tmpfreq = selectdata(tmpfreq{:}, 'param', 'fourierspctrm');
  %tmpfreq = combineplanar(cfg1,    tmpfreq);
  %tmpfreq = freqdescriptives(cfg2, tmpfreq);
  %cfg3.design = [ones(1,n1) ones(1,n2)*2];
  %tmpfreq.powspctrm = log10(tmpfreq.powspctrm);
  %tmpstat1 = freqstatistics(cfg3,   tmpfreq);
  %tmpstat1 = rmfield(tmpstat1, 'cfg');
  %tmpok    = double(cat(1, oktrials{mod(contrasts(cc,:)-1,5)+1}));
  %tmpok(tmpok==0) = nan;
  %tmpfreq.powspctrm = tmpfreq.powspctrm.*tmpok;
  %tmpstat2 = freqstatistics(cfg3,   tmpfreq);
  %tmpstat2 = rmfield(tmpstat2, 'cfg');
  cfg3.design = zeros(1,length(alldesign));
  cfg3.design(alldesign==contrasts(cc,1)) = 1;
  cfg3.design(alldesign==contrasts(cc,2)) = 2;
  tmpstat1 = freqstatistics(cfg3, freq);
  tmpstat1 = rmfield(tmpstat1, 'cfg');
  %tmpok    = double(cat(1, oktrials{[1:5 1:5]}));
  %tmpok(tmpok==0) = nan;
  tmpfreq  = freq;
  %tmpfreq.powspctrm = freq.powspctrm.*tmpok;
  n1 = sum(cfg3.design==1);
  n2 = sum(cfg3.design==2);
  n  = min(n1,n2);
  sel1 = randperm(n1);x1 = zeros(1,n1); x1(sel1(1:n)) = 1;
  sel2 = randperm(n2);x2 = zeros(1,n2); x2(sel2(1:n)) = 2;
  cfg3.design(cfg3.design==1) = x1;
  cfg3.design(cfg3.design==2) = x2;
  tmpstat2 = freqstatistics(cfg3, tmpfreq);
  tmpstat2 = rmfield(tmpstat2, 'cfg');

  stat(cc,1) = tmpstat1;
  stat(cc,2) = tmpstat2;
end
