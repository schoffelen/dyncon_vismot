global ft_default
ft_default.checksize = inf;

load vismot_parcels
load list

subjectname = list{s};
if s<3
  cd /project/3011085.03/analysis/mim/
else
  cd /project_ext/3010029/reproducescript/analysis/mim/
end


x = load(sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subjectname, 'previous'));
T = numel(x.trialnumber);

filename = sprintf('%s_mim_pre_all', subjectname);
mim_all = load([sprintf('%s/', subjectname), filename]);
c = mim_all.mim.mimspctrm;

filename = ['singletrial/', subjectname, '/', filename, '*'];
d = dir(filename);
dim = 6;

if dim~=6
  P3=1;
end
if dim~=16 && dim~=6
  P2=1;
end

freqs = 1:1:90;%0:0.5:119.5;
nfreq = 90;%239;

% Compute pseudo values and average within regions/frequencies of interest.
Ci = zeros(numel(d), dim, dim, nfreq);
trialnumber = zeros(T,1);
for k = 1:numel(d)
  k
  trialnumber(k,1) = str2num(char(extractBetween(d(k).name,"all_", ".mat")));
  m = load([d(k).folder,'/', d(k).name]);
  ci = m.mim.mimspctrm;
  Citmp = T*c-(T-1)*ci; % see Womelsdorf et al, 2007, Science
  for l=1:size(Citmp,3)
    Citmp2(:,:,l) = P3*P2*P*Citmp(:,:,l)*P'*P2'*P3';
  end
  if nfreq==4
    Ci(k,:,:,1)  = mean(Citmp2(:,:,nearest(freqs,8):nearest(freqs,12)),3);
    Ci(k,:,:,2)  = mean(Citmp2(:,:,nearest(freqs,12):nearest(freqs,30)),3);
    Ci(k,:,:,3)  = mean(Citmp2(:,:,nearest(freqs,30):nearest(freqs,50)),3);
    Ci(k,:,:,4)  = mean(Citmp2(:,:,nearest(freqs,50):nearest(freqs,70)),3);
    freqs = [10 21 40 60];
  elseif nfreq==90 || nfreq==80 || nfreq==120
    for g=1:nfreq
      Ci(k,:,:,g) = mean(Citmp2(:,:,2*g-1:2*g),3);
    end
  else
    Ci(k,:,:,:) = Citmp2(:,:,2:nfreq+1);
  end
  clear Citmp Citmp2 ci m
end

% find the condition number of every trial
trialinfo = zeros(size(x.trialinfo));
for k=1:numel(trialnumber)
  trialinfo(k,:) = x.trialinfo(x.trialinfo(:,2)==trialnumber(k),:);
end
conditions = trialinfo(:,end);

% filename = sprintf('/project/3011085.03/analysis/mim/singletrial/%s/%s_pseudomim', subjectname, subjectname);
% save(filename, 'Ci', 'conditions', 'trialnumber', 'trialinfo', 'nfreq', 'dim', 'freqs');

% remove outliers
[a1,a2,a3,a4] = size(Ci);
tmpdata.trial = reshape(Ci, [a1, 1, a2*a3*a4]);
tmpdata.dimord = 'rpt_chan_time';
tmpdata.label = {'chan01'};
tmpdata.time = 1:size(tmpdata.trial,3);
tmpdata.trialinfo = trialinfo;
cfg=[];
cfg.method = 'summary';
rej = ft_rejectvisual(cfg, tmpdata);

Ci = reshape(rej.trial, [size(rej.trial,1), a2, a3, a4]);
conditions = rej.trialinfo(:,end);
trialinfo = rej.trialinfo;

% remove neutral condition
x = find(conditions==5);
Ci(x,:,:,:)=[];
conditions(x)=[];
trialinfo(x,:) = [];

% equalize conditions
ntrials = min([sum(conditions==1), sum(conditions==2), sum(conditions==3), sum(conditions==4)]);
tmpidx = [];
for k=1:numel(unique(conditions))
  tmp = find(conditions==k);
  P = randperm(numel(tmp));
  tmpidx = [tmpidx; tmp(P(1:ntrials))];
end
Ci = Ci(tmpidx,:,:,:);
conditions = conditions(tmpidx);
trialinfo = trialinfo(tmpidx,:);

% Only take lower triangle of dim dimensions
sel = tril(ones(dim),-1)==1;
Citmp = reshape(Ci, [size(Ci,1), dim*dim, size(Ci,4)]);
Citmp = Citmp(:,sel,:);
% Citmp = reshape(Citmp, [size(Citmp,1), 1, size(Citmp,2)*size(Citmp,3)]);

data=[];
data.trial = Citmp;
data.time = 1:size(data.trial,3);
for k=1:size(data.trial,2)
  data.label{k} = sprintf('chan%0d', k);
end
% data.label = {'chan01'};
data.dimord = 'rpt_chan_time';

addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

cfg=[];
cfg.method = 'crossvalidate';
cfg.mva = {dml.standardizer dml.naive};%dml.one_against_rest('mva', dml.svm)}
cfg.statistic = {'confusion', 'accuracy'};
cfg.design = conditions;
cfg.type = 'nfold';
cfg.nfolds = 5;%numel(conditions);
cfg.resample = 0; % resamples conditions with fewer trials, and throws away oversampled conditions
stat = ft_timelockstatistics(cfg, data);
stat.statistic

%%
x(:,:,1) = mean(data.trial(:,:,8:12),3);
x(:,:,2) = mean(data.trial(:,:,12:30),3);
x(:,:,3) = mean(data.trial(:,:,30:50),3);
x(:,:,4) = mean(data.trial(:,:,50:70),3);
x = reshape(data.trial, numel(conditions), []);
rt = trialinfo(:,3);
mdl = fitglm(x, rt)

idx1 = conditions==1;
x1 = x(idx1);
rt1 = rt(idx1);
mdl1  = fitglm(x1, rt1)

idx2 = conditions==2;
x2 = x(idx2);
rt2 = rt(idx2);
mdl2  = fitglm(x2, rt2)

idx3 = conditions==3;
x3 = x(idx3);
rt3 = rt(idx3);
mdl3  = fitglm(x3, rt3)

idx4 = conditions==4;
x4 = x(idx4);
rt4 = rt(idx4);
mdl4  = fitglm(x4, rt4)

%%
% S(h) = stat.statistic.accuracy
figure; imagesc(stat.statistic.confusion)

numrandomization = 100;
for k=1:numrandomization
  dum{k} = conditions(randperm(numel(conditions)));
end

r=zeros(numrandomization,1);
for k=1:numrandomization
  k
  cfg.design = dum{k};
  randstat{k} = ft_timelockstatistics(cfg, data);
  r(k) = randstat{k}.statistic.accuracy;
end

filename = sprintf('/project/3011085.03/analysis/mim/singletrial/%s/%s_mimdecoding', subjectname, subjectname);
save(filename, 'data', 'conditions', 'cfg', 'stat', 'dim', 'nfreq', 'randstat', 'r');

if ~exist('splitlr', 'var'); splitlr=true; end
if splitlr
  data13 = data;
  x=find(conditions==1 | conditions==3);
  data13.trial=data.trial(x,:,:);
  conditions13 = conditions(x);
  cfg.design = conditions13;
  cfg.design(cfg.design==3)=2;
  stat13 = ft_timelockstatistics(cfg, data13);
  stat13.statistic
  
  data42 = data;
  x=find(conditions==4 | conditions==2);
  data42.trial=data.trial(x,:,:);
  conditions42 = conditions(x);
  cfg.design = conditions42;
  cfg.design(cfg.design==4)=1;
  stat42 = ft_timelockstatistics(cfg, data42);
  stat42.statistic
  
  % hemiflip condition 4 and 2, and pool with conditions 1 and 3.
  dataResp = data;
  x=find(conditions==4 | conditions==2);
  tmp42 = Ci(x,:,:,:);
  rev = [dim/2+1:dim 1:dim/2];
  tmp42 = tmp42(:, rev, rev,:);
  
  % Only take lower triangle of dim dimensions
  tmp42 = reshape(tmp42, [size(tmp42,1), dim*dim, size(tmp42,4)]);
  tmp42 = tmp42(:,sel,:);
  
  dataResp.trial(x,:,:) = tmp42;
  cfg.design = conditions;
  cfg.design(cfg.design==4)=1;
  cfg.design(cfg.design==3)=2;
  statResp = ft_timelockstatistics(cfg, dataResp);
  statResp.statistic
  
  save(filename, 'stat13', 'stat42', 'statResp', '-append');
  figure; imagesc(stat13.statistic.confusion)
  figure; imagesc(stat42.statistic.confusion)
  figure; imagesc(statResp.statistic.confusion)
end



