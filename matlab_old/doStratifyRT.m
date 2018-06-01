function [stat13, stat42, stat13b, stat42b] = doStratifyRT(subject, flag)

if nargin<2,
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
if flag==0,
  load([subject.name,'mtmfft004bnew']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned004bnew']);
elseif flag==2,
  load([subject.name,'mtmfft_aligned004bcnew']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
rtall = [];
for k = 1:5
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+5) = length(freqpst{k}.cumtapcnt);
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
oklabel = gradlabel(ok==5);
for k = 1:5
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', 10);
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', 10);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end
freq = selectdata(freqpre{:}, freqpst{:}, 'param', 'fourierspctrm');
clear freqpre freqpst

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);
%1:5 are the baselines for conditions 1:5
%6:10 are the post stim intervals for conditions 1:5

%congruency and RT glm for response left
%activation condition 1 vs 3

%stratify for RT
sel1   = btrl(1):etrl(1);
sel2   = btrl(3):etrl(3);
rt1    = freq.rt(sel1);
rt2    = freq.rt(sel2);
output = stratify([],rt1(:)',rt2(:)');
sel1   = sel1(find(isfinite(output{1})));
sel2   = sel2(find(isfinite(output{2})));


%congruency and RT glm for response right
%activation condition 4 vs 2

%stratify for RT
sel1   = btrl(4):etrl(4);
sel2   = btrl(2):etrl(2);
rt1    = freq.rt(sel1);
rt2    = freq.rt(sel2);
output = stratify([],rt1(:)',rt2(:)');
sel1   = sel1(find(isfinite(output{1})));
sel2   = sel2(find(isfinite(output{2})));

