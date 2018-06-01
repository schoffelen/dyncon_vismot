function [stat13, stat42, stat13b, stat42b] = doFreqanalysisPrePstStratified(subject, frequency, flag)

%Try to do the something similar as doSourceanalysisDICSprepst3
%
%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's baseline and post
%stimulus interval are matched, and number of samples across conditions which are contrasted are matched)
%spatial filters are computed on all trials concatenated
%output:
%GLM analysis
%based on data stratified for RT

%stat13 stat42: act congruency effect in field stat and stat2
%stat13 stat42: act rt effect in field statrt stat2rt
%stat13b stat42b: as in act but then baselines

if nargin<3,
  flag = 0;
end

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

if length(subject.runnames)>1,
  for k = 1:length(subject.runnames)
    fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{k},subject.datafile];
    event{k}.event = read_event(fname);
  end
else
  event = read_event(fname);
end

cd([subject.pathname,'freq/']);
freqx = 10;
if flag==0 && freqx<37,
  load([subject.name,'mtmfft004bnew']);
elseif flag==0,
  load([subject.name,'mtmfft004bnew']);
  for k = 1:5
    trl1{k,1} = findcfg(freqpre{k}.cfg, 'trl');
    trl2{k,1} = findcfg(freqpst{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl1{k};
    freqpst{k}.cfg.trl = trl2{k};
  end
elseif flag==1 && freqx<37,
  load([subject.name,'mtmfft_aligned004bnew']);
elseif flag==2 && freqx<37,
  load([subject.name,'mtmfft_aligned004bcnew']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned004bnew']);
  for k = 1:5
    trl1{k,1} = findcfg(freqpre{k}.cfg, 'trl');
    trl2{k,1} = findcfg(freqpst{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft_aligned012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl1{k};
    freqpst{k}.cfg.trl = trl2{k};
  end
elseif flag==2,
  load([subject.name,'mtmfft_aligned012bcnew']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
rtall = [];
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+4) = length(freqpst{k}.cumtapcnt);
  if k==5.
    rt = nan+zeros(length(freqpre{k}.cumtapcnt),1);
    freqpre{k}.rt = rt;
    rt = nan+zeros(length(freqpst{k}.cumtapcnt),1);
    freqpst{k}.rt = rt;
  else
    rt = data2rt(freqpre{k}, event, 2);
    freqpre{k}.rt = rt./fsample;
    freqpst{k}.rt = rt./fsample;
  end
  rtall = [rtall;rt./fsample]; %rt for condition 5 does not mean anything
end
oklabel = gradlabel(ok==4);
for k = 1:4
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', frequency);
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

load([subject.pathname,'stratifyRT/',subject.name,'stratifyRT.mat']);
for k = 1:4
  freqpre{k} = selectdata(freqpre{k}, 'rpt', find(isfinite(output{k})));
  freqpst{k} = selectdata(freqpst{k}, 'rpt', find(isfinite(output{k})));
  freqpre{k}.rt = freqpre{k}.rt(find(isfinite(output{k})));
  freqpst{k}.rt = freqpst{k}.rt(find(isfinite(output{k})));
  ntrl(k) = sum(isfinite(output{k}));
  ntrl(k+4) = sum(isfinite(output{k}));
end
freq{1} = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
freq{2} = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
freq    = selectdata(freq{:}, 'param', 'fourierspctrm');
rtall   = cat(2,output{:});
rtall   = rtall(isfinite(rtall));
clear freqpre freqpst

cfgp  = [];
cfgp.planarmethod  = 'sincos';
cfgp.neighbourdist = 0.037;
freq  = megplanar(cfgp, freq);
freq = checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
cfgp2 = [];
cfgp2.combinemethod = 'sum';
freq  = combineplanar(cfgp2, freq);

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);
%1:4 are the baselines for conditions 1:4
%5:8 are the post stim intervals for conditions 1:4

freq.rt = log([rtall(:);rtall(:)]);

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.parameter = 'powspctrm';
cfg2.statistic = 'glm';
cfg2.glm.statistic = 'T';

%log transform powerspctrm
freq.powspctrm = log(freq.powspctrm);

%congruency and RT glm for response left
%activation condition 1 vs 3

%stratify for RT
sel1   = btrl(5):etrl(5);
sel2   = btrl(7):etrl(7);

tmpsd = selectdata(freq, 'rpt', [sel1 sel2]);
tmpsd.powspctrm = standardise(tmpsd.powspctrm, 1);
n1    = numel(sel1);
n2    = numel(sel2);
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
%%orthogonalise with respect to the condition
%rt1   = rt1(:)' - mean(rt1(:));
%rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(3,:) = [freq.cumtapcnt([sel1 sel2])];
cfg2.design(2,:) = ([rt1(:)' rt2(:)']);
cfg2.design      = orthogonalise(cfg2.design')';

cfg2.glm.contrast  = [1 0 0];
stat13     = freqstatistics(cfg2, tmpsd);
stat13.cfg   = [];
cfg2.glm.contrast = [0 0 1]; %correlation with rt
stat13rt       = freqstatistics(cfg2, tmpsd);
stat13.stat2   = stat13rt.stat;
cfg2.glm.contrast  =[0 1 0];
stat13nt       = freqstatistics(cfg2, tmpsd);
stat13.stat3   = stat13nt.stat;

%congruency and RT glm for response right
%activation condition 4 vs 2

%stratify for RT
sel1   = btrl(8):etrl(8);
sel2   = btrl(6):etrl(6);

tmpsd = selectdata(freq, 'rpt', [sel1 sel2]);
tmpsd.powspctrm = standardise(tmpsd.powspctrm, 1);
n1    = numel(sel1);
n2    = numel(sel2);
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
%%orthogonalise with respect to the condition
%rt1   = rt1(:)' - mean(rt1(:));
%rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(3,:) = [freq.cumtapcnt([sel1 sel2])];
cfg2.design(2,:) = ([rt1(:)' rt2(:)']);
cfg2.design      = orthogonalise(cfg2.design')';
cfg2.glm.contrast  = [1 0 0];
stat42     = freqstatistics(cfg2, tmpsd);
stat42.cfg   = [];
cfg2.glm.contrast = [0 0 1]; %correlation with rt
stat42rt       = freqstatistics(cfg2, tmpsd);
stat42.stat2   = stat42rt.stat;
cfg2.glm.contrast  =[0 1 0];
stat42nt       = freqstatistics(cfg2, tmpsd);
stat42.stat3   = stat42nt.stat;

%congruency and RT glm for response left baseline
%baseline condition 1 vs 3

%stratify for RT
sel1   = btrl(1):etrl(1);
sel2   = btrl(3):etrl(3);

tmpsd = selectdata(freq, 'rpt', [sel1 sel2]);
tmpsd.powspctrm = standardise(tmpsd.powspctrm, 1);
n1    = numel(sel1);
n2    = numel(sel2);
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
%%orthogonalise with respect to the condition
%rt1   = rt1(:)' - mean(rt1(:));
%rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(3,:) = [freq.cumtapcnt([sel1 sel2])];
cfg2.design(2,:) = ([rt1(:)' rt2(:)']);
cfg2.design      = orthogonalise(cfg2.design')';
cfg2.glm.contrast  = [1 0 0];
stat13b     = freqstatistics(cfg2, tmpsd);
stat13b.cfg   = [];
cfg2.glm.contrast = [0 0 1]; %correlation with rt
stat13brt       = freqstatistics(cfg2, tmpsd);
stat13b.stat2   = stat13brt.stat;
cfg2.glm.contrast  =[0 1 0];
stat13bnt       = freqstatistics(cfg2, tmpsd);
stat13b.stat3   = stat13bnt.stat;

%congruency and RT glm for response right baseline
%baseline condition 4 vs 2

%stratify for RT
sel1   = btrl(4):etrl(4);
sel2   = btrl(2):etrl(2);

tmpsd = selectdata(freq, 'rpt', [sel1 sel2]);
tmpsd.powspctrm = standardise(tmpsd.powspctrm, 1);
n1    = numel(sel1);
n2    = numel(sel2);
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
%%orthogonalise with respect to the condition
%rt1   = rt1(:)' - mean(rt1(:));
%rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(3,:) = [freq.cumtapcnt([sel1 sel2])];
cfg2.design(2,:) = ([rt1(:)' rt2(:)']);
cfg2.design      = orthogonalise(cfg2.design')';
cfg2.glm.contrast  = [1 0 0];
stat42b     = freqstatistics(cfg2, tmpsd);
stat42b.cfg   = [];
cfg2.glm.contrast = [0 0 1]; %correlation with rt
stat42brt       = freqstatistics(cfg2, tmpsd);
stat42b.stat2   = stat42brt.stat;
cfg2.glm.contrast  =[0 1 0];
stat42bnt       = freqstatistics(cfg2, tmpsd);
stat42b.stat3   = stat42bnt.stat;


function [x] = orthogonalise(input)

input = standardise(input,1);
x     = input(:,1);

for i = 2:size(input,2)
  y   = input(:,i);
  y   = y-x*(pinv(x'*x)*x'*y);
  if any(y), x = [x y]; end
end
