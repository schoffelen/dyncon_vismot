function [output, input, binaxis] = stratifyRT(subject,  flag)

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

frequency = 10;
cd([subject.pathname,'freq/']);
if flag==0 && frequency<37,
  load([subject.name,'mtmfft004bnew']);
elseif flag==0,
  load([subject.name,'mtmfft004bnew']);
  for k = 1:5
    trl{k,1} = findcfg(freqpre{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl{k};
  end
elseif flag==1 && frequency<37,
  load([subject.name,'mtmfft_aligned004bnew']);
elseif flag==2 && frequency<37,
  load([subject.name,'mtmfft_aligned004bcnew']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned004bnew']);
  for k = 1:5
    trl{k,1} = findcfg(freqpre{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft_aligned012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl{k};
  end
elseif flag==2,
  load([subject.name,'mtmfft_aligned012bcnew']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
rtall = [];
for k = 1:5
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpst{k}.cumtapcnt);
  if k==5.
    rt = nan+zeros(length(freqpst{k}.cumtapcnt),1);
    freqpst{k}.rt = rt;
  else
    rt = data2rt(freqpre{k}, event, 2);
    freqpst{k}.rt = rt./fsample;
  end
  rtall = [rtall;rt./fsample]; %rt for condition 5 does not mean anything
end

rt1    = freqpst{1}.rt';
rt2    = freqpst{2}.rt';
rt3    = freqpst{3}.rt';
rt4    = freqpst{4}.rt';
cfg.equalbinavg = 'no';
cfg.numbin = 8;
[output, binaxis] = stratify(cfg,rt1,rt2,rt3,rt4);
input{1,1} = rt1;
input{1,2} = rt2;
input{1,3} = rt3;
input{1,4} = rt4;
