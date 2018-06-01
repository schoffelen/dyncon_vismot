function [freqpre,freqpst] = doFreqanalysisBalancePrePostResplock(subject,smoothing,foilim,flag);

%this function does something similar to doFreqanalysisBalancePrePost2 but now response locked
%in addition to matching the number of trials across the relevant conditions, and the number of samples approximately within each condition (baseline versus stimulation), now also do a paired matching (baseline versus stimulation) for the number of samples, and keep track of rt

if smoothing==0,
  taper = 'hanning';
  smoothing = 4;
else
  taper = 'dpss';
end

if nargin<4,
  flag = 0;
end

cd(subject.pathname);
cd('data');
if flag==0,
  load([subject.name,'data']);
elseif flag==1,
  load([subject.name,'data_aligned']);
elseif flag==2,
  load([subject.name,'data_alignedclean']);
end

for k = 1:5
  warning off;
  if k==1,
    data = data1;
  elseif k==2,
    data = data2;
  elseif k==3,
    data = data3;
  elseif k==4,
    data = data4;
  elseif k==5,
    data = data5;
  end
  data              = struct2double(data);
  if flag==2, data  = fixResample(data); end

  cfg           = [];
  % this is the original
  %cfg.toilim    = [0.2 0.7-1./data.fsample];
  % this is the new one
  cfg.toilim    = [0.2 1.5];
  cfg.minlength = 0.25;
  datapst       = redefinetrial(cfg, data);
  cfg.toilim    = [-1.5 0-1/256];
  datapre       = redefinetrial(cfg, data);
  clear data;

  %match the pre and pst trials so that data can be stratified based on 
  %baseline power
  trlold  = datapre.cfg.trlold;
  trl     = datapre.cfg.trl;
  origvec = zeros(size(trlold,1),2);
  runold  = [0;find(diff(trlold(:,1))<0)]+1;
  runold(:,2) = [runold(2:end,1)-1;size(trlold,1)];
  run     = [0;find(diff(trl(:,1))<0)]+1;
  run(:,2) = [run(2:end,1)-1;size(trl,1)];
  for kk = 1:size(trl,1)
    block = find(kk<=run(:,2),1,'first');
    ii    = runold(block,1):runold(block,2);
    indx  = find(trlold(ii,1)<=trl(kk,1) & trlold(ii,2)>=trl(kk,2)) + runold(block,1) - 1;
    origvec(indx,1) = kk;
  end
  trlold  = datapst.cfg.trlold;
  trl     = datapst.cfg.trl;
  runold  = [0;find(diff(trlold(:,1))<0)]+1;
  runold(:,2) = [runold(2:end,1)-1;size(trlold,1)];
  run     = [0;find(diff(trl(:,1))<0)]+1;
  run(:,2) = [run(2:end,1)-1;size(trl,1)];
  for kk = 1:size(trl,1)
    block = find(kk<=run(:,2),1,'first');
    ii    = runold(block,1):runold(block,2);
    indx  = find(trlold(ii,1)<=trl(kk,1) & trlold(ii,2)>=trl(kk,2)) + runold(block,1) - 1;
    origvec(indx,2) = kk;
  end
  ok = origvec(:,1) > 0 & origvec(:,2) > 0;
  origvec = origvec(ok,:);

  cfg         = [];
  cfg.detrend = 'yes';
  cfg.trials  = origvec(:,1);
  datapre     = preprocessing(cfg, datapre);
  cfg.trials  = origvec(:,2);
  datapst     = preprocessing(cfg, datapst);
  Npre(k,1) = length(datapre.trial);
  Npst(k,1) = length(datapst.trial);
  pre{k}    = datapre;
  pst{k}    = datapst;
  clear datapre datapst;
end

%match number of trials and length
%match 1 with 3 and 4 with 2, 5 is ok as it is

smp1    = cellfun('size',pst{1}.trial,2);
smp2    = cellfun('size',pst{3}.trial,2);
for k = 1:length(smp1)
  M = min(smp1(k),128);
  pst{1}.trial{k} = pst{1}.trial{k}(:,end+1-M:end);
  pst{1}.time{k}  = pst{1}.time{k}(end+1-M:end);
end
for k = 1:length(smp2)
  M = min(smp2(k),128);
  pst{3}.trial{k} = pst{3}.trial{k}(:,end+1-M:end);
  pst{3}.time{k}  = pst{3}.time{k}(end+1-M:end);
end
smp1    = cellfun('size',pst{1}.trial,2);
smp2    = cellfun('size',pst{3}.trial,2);
[indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2);
for k = 1:length(indx1)
  %this is original stim locked
  %tmptrial1{1,k} = blc(pst{1}.trial{indx1(k)}(:,1:nsmp1(k)));
  %tmptrial2{1,k} = blc(pst{3}.trial{indx2(k)}(:,1:nsmp2(k)));
  %tmptime1{1,k}  = (pst{1}.time{indx1(k)}(1:nsmp1(k)));
  %tmptime2{1,k}  = (pst{3}.time{indx2(k)}(1:nsmp2(k)));
  
  %this is resp locked
  tmptrial1{1,k} = blc(pst{1}.trial{indx1(k)}(:,end+1-nsmp1(k):end));
  tmptrial2{1,k} = blc(pst{3}.trial{indx2(k)}(:,end+1-nsmp2(k):end));
  tmptime1{1,k}  = (pst{1}.time{indx1(k)}(end+1-nsmp1(k):end));
  tmptime2{1,k}  = (pst{3}.time{indx2(k)}(end+1-nsmp2(k):end));
  b1 = min([size(pre{1}.trial{indx1(k)},2),nsmp1(k)]);
  b2 = min([size(pre{3}.trial{indx2(k)},2),nsmp2(k)]);
  
  %this is stim locked
  %tmptrial3{1,k} = blc(pre{1}.trial{indx1(k)}(:,1:b1));
  %tmptrial4{1,k} = blc(pre{3}.trial{indx2(k)}(:,1:b2));
  %tmptime3{1,k}  = (pre{1}.time{indx1(k)}(1:b1));
  %tmptime4{1,k}  = (pre{3}.time{indx2(k)}(1:b2));
  
  %this is resp locked
  tmptrial3{1,k} = blc(pre{1}.trial{indx1(k)}(:,end+1-b1:end));
  tmptrial4{1,k} = blc(pre{3}.trial{indx2(k)}(:,end+1-b2:end));
  tmptime3{1,k}  = (pre{1}.time{indx1(k)}(end+1-b1:end));
  tmptime4{1,k}  = (pre{3}.time{indx2(k)}(end+1-b2:end));
end
pst{1}.trial = tmptrial1;
pst{1}.time  = tmptime1;
trl = findcfg(pst{1}.cfg,'trl');pst{1}.cfg.trl = trl(indx1,:);
pst{3}.trial = tmptrial2;
pst{3}.time  = tmptime2;
trl = findcfg(pst{3}.cfg,'trl');pst{3}.cfg.trl = trl(indx2,:);
clear tmptrial1 tmptime1 tmptrial2 tmptime2;
pre{1}.trial = tmptrial3;
pre{1}.time  = tmptime3;
trl = findcfg(pre{1}.cfg,'trl');pre{1}.cfg.trl = trl(indx1,:);
pre{3}.trial = tmptrial4;
pre{3}.time  = tmptime4;
trl = findcfg(pre{3}.cfg,'trl');pre{3}.cfg.trl = trl(indx2,:);
clear tmptrial3 tmptime3 tmptrial4 tmptime4;

smp1    = cellfun('size',pst{2}.trial,2);
smp2    = cellfun('size',pst{4}.trial,2);
for k = 1:length(smp1)
  M = min(smp1(k),128);
  pst{2}.trial{k} = pst{2}.trial{k}(:,end+1-M:end);
  pst{2}.time{k}  = pst{2}.time{k}(end+1-M:end);
end
for k = 1:length(smp2)
  M = min(smp2(k),128);
  pst{4}.trial{k} = pst{4}.trial{k}(:,end+1-M:end);
  pst{4}.time{k}  = pst{4}.time{k}(end+1-M:end);
end
smp1    = cellfun('size',pst{2}.trial,2);
smp2    = cellfun('size',pst{4}.trial,2);
[indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2);
for k = 1:length(indx1)
  %this is stim locked
  %tmptrial1{1,k} = blc(pst{2}.trial{indx1(k)}(:,1:nsmp1(k)));
  %tmptrial2{1,k} = blc(pst{4}.trial{indx2(k)}(:,1:nsmp2(k)));
  %tmptime1{1,k}  = (pst{2}.time{indx1(k)}(1:nsmp1(k)));
  %tmptime2{1,k}  = (pst{4}.time{indx2(k)}(1:nsmp2(k)));
  
  %this is resp locked
  tmptrial1{1,k} = blc(pst{2}.trial{indx1(k)}(:,end+1-nsmp1(k):end));
  tmptrial2{1,k} = blc(pst{4}.trial{indx2(k)}(:,end+1-nsmp2(k):end));
  tmptime1{1,k}  = (pst{2}.time{indx1(k)}(end+1-nsmp1(k):end));
  tmptime2{1,k}  = (pst{4}.time{indx2(k)}(end+1-nsmp2(k):end));
  b1 = min([size(pre{2}.trial{indx1(k)},2),nsmp1(k)]);
  b2 = min([size(pre{4}.trial{indx2(k)},2),nsmp2(k)]);
 
  %this is stim locked
  %tmptrial3{1,k} = blc(pre{2}.trial{indx1(k)}(:,1:b1));
  %tmptrial4{1,k} = blc(pre{4}.trial{indx2(k)}(:,1:b2));
  %tmptime3{1,k}  = (pre{2}.time{indx1(k)}(1:b1));
  %tmptime4{1,k}  = (pre{4}.time{indx2(k)}(1:b2));
  
  %this is resp locked
  tmptrial3{1,k} = blc(pre{2}.trial{indx1(k)}(:,end+1-b1:end));
  tmptrial4{1,k} = blc(pre{4}.trial{indx2(k)}(:,end+1-b2:end));
  tmptime3{1,k}  = (pre{2}.time{indx1(k)}(end+1-b1:end));
  tmptime4{1,k}  = (pre{4}.time{indx2(k)}(end+1-b2:end));
end
pst{2}.trial = tmptrial1;
pst{2}.time  = tmptime1;
trl = findcfg(pst{2}.cfg,'trl');pst{2}.cfg.trl = trl(indx1,:);
pst{4}.trial = tmptrial2;
pst{4}.time  = tmptime2;
trl = findcfg(pst{4}.cfg,'trl');pst{4}.cfg.trl = trl(indx2,:);
clear tmptrial1 tmptime1 tmptrial2 tmptime2;
pre{2}.trial = tmptrial3;
pre{2}.time  = tmptime3;
trl = findcfg(pre{2}.cfg,'trl');pre{2}.cfg.trl = trl(indx1,:);
pre{4}.trial = tmptrial4;
pre{4}.time  = tmptime4;
trl = findcfg(pre{4}.cfg,'trl');pre{4}.cfg.trl = trl(indx2,:);
clear tmptrial3 tmptime3 tmptrial4 tmptime4;

smp = cellfun('size',pst{5}.trial,2);
for k = 1:length(smp)
  nsmp = min([smp(k),size(pre{5}.trial{k},2)]);
  pre{5}.trial{k} = blc(pre{5}.trial{k}(:,1:nsmp));
  pre{5}.time{k}  = pre{5}.time{k}(1:nsmp);
end

cfg         = [];
cfg.method  = 'mtmfft';
cfg.output  = 'fourier';
cfg.pad     = 512./pre{1}.fsample; %explicitly make nfft 256
cfg.foilim  = foilim;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
cfg.channel = 'MEG';
for k =1:4
  freqpre{k} = selectdata(freqanalysis(cfg, pre{k}),'foilim',foilim(1):2:foilim(2));
  freqpst{k} = selectdata(freqanalysis(cfg, pst{k}),'foilim',foilim(1):2:foilim(2));
end

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
