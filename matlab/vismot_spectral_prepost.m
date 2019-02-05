function [freqpre,freqpst] = vismot_spectral_prepost(subject,varargin)

% This function is based on doFreqanalysisMtmfft in the matlab_old folder
%
% do mtmfft analysis on data stratified for RT and balanced number of samples
% stimulus locked

smoothing = ft_getopt(varargin, 'smoothing', 4);
foilim    = ft_getopt(varargin, 'foilim', [0 80]);
output    = ft_getopt(varargin, 'output', 'pow');
doplanar  = istrue(ft_getopt(varargin, 'doplanar', false));

if smoothing==0
  taper = 'hanning';
else
  taper = 'dpss';
end

% load the data
load(fullfile(subject.pathname,'data',[subject.name,'data']));

pre = cell(1,5);
pst = cell(1,5);
for k = 1:5
  if k==1
    data = data1;
  elseif k==2
    data = data2;
  elseif k==3
    data = data3;
  elseif k==4
    data = data4;
  elseif k==5
    data = data5;
  end
  
  cfg           = [];
  cfg.toilim    = [0.2 0.7-1./data.fsample];
  cfg.minlength = 0.25;
  datapst       = ft_redefinetrial(cfg, data);
  cfg.toilim    = [-0.5 0-1./data.fsample];
  datapre       = ft_redefinetrial(cfg, data);
  clear data;
  
  cfg         = [];
  cfg.detrend = 'yes';
  datapre     = ft_preprocessing(cfg, datapre);
  datapst     = ft_preprocessing(cfg, datapst);
  pre{k}    = datapre;
  pst{k}    = datapst;
  clear datapre datapst;
end

%match number of trials and length based on the post data segments
%match 1 with 3 and 4 with 2, 5 is ok as it is

smp1    = cellfun('size',pst{1}.trial,2);
smp2    = cellfun('size',pst{3}.trial,2);
[indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2); % this function is in private

tmptrial1 = cell(1,numel(indx1));
tmptrial2 = cell(1,numel(indx1));
tmptrial3 = cell(1,numel(indx1));
tmptrial4 = cell(1,numel(indx1));
tmptime1 = cell(1,numel(indx1));
tmptime2 = cell(1,numel(indx1));
tmptime3 = cell(1,numel(indx1));
tmptime4 = cell(1,numel(indx1));

sel1 = false(numel(indx1),1);
sel2 = false(numel(indx2),1);
for k = 1:length(indx1)
  tmptrial1{1,k} = ft_preproc_baselinecorrect(pst{1}.trial{indx1(k)}(:,1:nsmp1(k)));
  tmptrial2{1,k} = ft_preproc_baselinecorrect(pst{3}.trial{indx2(k)}(:,1:nsmp2(k)));
  tmptime1{1,k}  = (pst{1}.time{indx1(k)}(1:nsmp1(k)));
  tmptime2{1,k}  = (pst{3}.time{indx2(k)}(1:nsmp2(k)));
  
  try
    b1             = min([size(pre{1}.trial{indx1(k)},2),nsmp1(k)]);
    tmptrial3{1,k} = ft_preproc_baselinecorrect(pre{1}.trial{indx1(k)}(:,1:b1));
    tmptime3{1,k}  = (pre{1}.time{indx1(k)}(1:b1));
    sel1(k) = true;
  end
  try
    b2             = min([size(pre{3}.trial{indx2(k)},2),nsmp2(k)]);
    tmptrial4{1,k} = ft_preproc_baselinecorrect(pre{3}.trial{indx2(k)}(:,1:b2));
    tmptime4{1,k}  = (pre{3}.time{indx2(k)}(1:b2));
    sel2(k) = true;
  end
end
pst{1}.trial = tmptrial1;
pst{1}.time  = tmptime1;
pst{1}.trialinfo = pst{1}.trialinfo(indx1,:);
pst{3}.trial = tmptrial2;
pst{3}.time  = tmptime2;
pst{3}.trialinfo = pst{3}.trialinfo(indx2,:);
clear tmptrial1 tmptime1 tmptrial2 tmptime2;
pre{1}.trial = tmptrial3(sel1);
pre{1}.time  = tmptime3(sel1);
pre{1}.trialinfo = pre{1}.trialinfo(indx1(sel1),:);
pre{3}.trial = tmptrial4(sel2);
pre{3}.time  = tmptime4(sel2);
pre{3}.trialinfo = pre{3}.trialinfo(indx2(sel2),:);
clear tmptrial3 tmptime3 tmptrial4 tmptime4;

smp1    = cellfun('size',pst{2}.trial,2);
smp2    = cellfun('size',pst{4}.trial,2);
[indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2);

tmptrial1 = cell(1,numel(indx1));
tmptrial2 = cell(1,numel(indx1));
tmptrial3 = cell(1,numel(indx1));
tmptrial4 = cell(1,numel(indx1));
tmptime1 = cell(1,numel(indx1));
tmptime2 = cell(1,numel(indx1));
tmptime3 = cell(1,numel(indx1));
tmptime4 = cell(1,numel(indx1));

sel1 = false(numel(indx1),1);
sel2 = false(numel(indx2),1);
for k = 1:length(indx1)
  tmptrial1{1,k} = ft_preproc_baselinecorrect(pst{2}.trial{indx1(k)}(:,1:nsmp1(k)));
  tmptrial2{1,k} = ft_preproc_baselinecorrect(pst{4}.trial{indx2(k)}(:,1:nsmp2(k)));
  tmptime1{1,k}  = (pst{2}.time{indx1(k)}(1:nsmp1(k)));
  tmptime2{1,k}  = (pst{4}.time{indx2(k)}(1:nsmp2(k)));
  
  try
    b1 = min([size(pre{2}.trial{indx1(k)},2),nsmp1(k)]);
    tmptrial3{1,k} = ft_preproc_baselinecorrect(pre{2}.trial{indx1(k)}(:,1:b1));
    tmptime3{1,k}  = (pre{2}.time{indx1(k)}(1:b1));
    sel1(k) = true;
  end
  try
    b2 = min([size(pre{4}.trial{indx2(k)},2),nsmp2(k)]);
    tmptrial4{1,k} = ft_preproc_baselinecorrect(pre{4}.trial{indx2(k)}(:,1:b2));
    tmptime4{1,k}  = (pre{4}.time{indx2(k)}(1:b2));
    sel2(k) = true;
  end
  
end
pst{2}.trial = tmptrial1;
pst{2}.time  = tmptime1;
pst{2}.trialinfo = pst{2}.trialinfo(indx1,:);
pst{4}.trial = tmptrial2;
pst{4}.time  = tmptime2;
pst{4}.trialinfo = pst{4}.trialinfo(indx2,:);
clear tmptrial1 tmptime1 tmptrial2 tmptime2;
pre{2}.trial = tmptrial3(sel1);
pre{2}.time  = tmptime3(sel1);
pre{2}.trialinfo = pre{2}.trialinfo(indx1(sel1),:);
pre{4}.trial = tmptrial4(sel2);
pre{4}.time  = tmptime4(sel2);
pre{4}.trialinfo = pre{4}.trialinfo(indx2(sel2),:);
clear tmptrial3 tmptime3 tmptrial4 tmptime4;

cfg         = [];
cfg.method  = 'mtmfft';
cfg.output  = output;
cfg.pad     = 300./pre{1}.fsample; %explicitly make nfft 300
cfg.foilim  = foilim;
cfg.taper   = taper;
cfg.tapsmofrq = smoothing;
cfg.channel   = pre{1}.label;
cfg.keeptrials = 'yes';

% convert to synthetic planar gradient
if doplanar
  cfgn = [];
  cfgn.method   = 'template';
  cfgn.template = 'bti248_neighb.mat';
  neighbours = ft_prepare_neighbours(cfgn);
  
  cfgp = [];
  cfgp.method = 'sincos';
  cfgp.neighbours = neighbours;
  for k = 1:5
    pre{k} = ft_megplanar(cfgp, pre{k});
    pst{k} = ft_megplanar(cfgp, pst{k});
  end
end

for k = 1:5
  freqpre(k) = ft_freqanalysis(rmfield(cfg, 'channel'), pre{k});
  freqpst(k) = ft_freqanalysis(rmfield(cfg, 'channel'), pst{k});
end

% combine planar gradients
if doplanar
  for k = 1:5
    freqpre(k) = ft_combineplanar([], freqpre(k));
    freqpst(k) = ft_combineplanar([], freqpst(k));
  end
end

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
