load vismot_parcels
cd /project/3011085.03/analysis/mim
load list

subjectname = list{2};

x = load(sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subjectname, 'previous'));
T = numel(x.trialnumber);

filename = sprintf('%s_mim_pre_all', subjectname);
mim_all = load(filename);
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

freqs = 0.5:0.5:120;
nfreq = 80;

% Compute pseudo values and average within regions/frequencies of interest.
Ci = zeros(numel(d), dim, dim, nfreq);
for k = 1:numel(d)
  k
  trialnumber(k) = str2num(char(extractBetween(d(k).name,"all_", ".mat")));
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
  elseif nfreq==80
    Ci(k,:,:,:) = Citmp2(:,:,3:2:2*(nfreq+1));
    freqs = 1:1:80;
  else
    Ci(k,:,:,:) = Citmp2(:,:,2:nfreq+1);
  end
    clear Citmp Citmp2 ci m
end

% find the condition number of every trial
[initnumber,idx] = unique(x.trialnumber);
tmptrialinfo = x.trialinfo(idx,:);
conditions = zeros(numel(idx),1);
trialinfo = zeros(size(tmptrialinfo));
for k=1:numel(idx)
  conditions(k,1) = tmptrialinfo((initnumber==trialnumber(k)), end);
  trialinfo(k,:) = tmptrialinfo((initnumber==trialnumber(k)), :);
end

% filename = sprintf('/project/3011085.03/analysis/mim/singletrial/%s/%s_pseudomim', subjectname, subjectname);
% save(filename, 'Ci', 'conditions', 'trialnumber', 'trialinfo', 'nfreq', 'dim', 'freqs');

% remove outliers
[a1,a2,a3,a4] = size(Ci);
tmp = reshape(Ci, [a1, a2*a3*a4]);
tmpdata.trial = reshape(tmp, [size(tmp,1), 1, size(tmp,2)]);
tmpdata.dimord = 'rpt_chan_time';
tmpdata.label = {'chan01'};
tmpdata.time = 1:size(tmpdata.trial,3);
tmpdata.trialinfo = conditions;
cfg=[];
cfg.method = 'summary';
rej = ft_rejectvisual(cfg, tmpdata);

Ci = reshape(rej.trial, [size(rej.trial,1), a2, a3, a4]); 
conditions = rej.trialinfo;

% remove neutral condition
x = find(conditions==5);
Ci(x,:,:,:)=[];
conditions(x)=[];

% equalize conditions
ntrials = min([sum(conditions==1), sum(conditions==2), sum(conditions==3), sum(conditions==4)]);
tmpidx = [];
for k=1:4
  tmp = find(conditions==k);
  P = randperm(sum(conditions==k));
  tmpidx = [tmpidx; tmp(P(1:ntrials))];
end
Ci = Ci(tmpidx,:,:,:);
conditions = conditions(tmpidx);

% Only take lower triangle of dim dimensions
sel = tril(ones(dim),-1)==1;
Citmp = reshape(Ci, [size(Ci,1), dim*dim, size(Ci,4)]);
Citmp = Citmp(:,sel,:);
Citmp = reshape(Citmp, [size(Citmp,1), 1, size(Citmp,2)*size(Citmp,3)]);

data=[];
data.trial = Citmp;
data.time = 1:size(data.trial,3);
data.label = {'chan01'};
data.dimord = 'rpt_chan_time';

addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

cfg=[];
cfg.method = 'crossvalidate';
cfg.mva = {dml.standardizer dml.naive};
cfg.statistic = {'confusion', 'accuracy'};
cfg.design = conditions;
cfg.type = 'nfold';
cfg.nfolds = 5;%numel(conditions);
cfg.resample = 0; % resamples conditions with fewer trials, and throws away oversampled conditions
stat = ft_timelockstatistics(cfg, data);
stat.statistic


% filename = sprintf('/project/3011085.03/analysis/mim/singletrial/%s/%s_mimdecoding', subjectname, subjectname);
% save(filename, 'data', 'conditions', 'cfg', 'stat', 'dim', 'nfreq');

%%
numrandomization = 100;
r=zeros(numrandomization,1);
for k=1:numrandomization
  k
  cfg.design = conditions(randperm(numel(conditions)));
  randstat{k} = ft_timelockstatistics(cfg, data);
  r(k) = randstat{k}.statistic.accuracy;
end
%}
%{
data13 = data;
x=find(conditions==1 | conditions==3);
data13.trial=data.trial(x,:,:);
conditions13 = conditions(x);
cfg.design = conditions13;
cfg.design(cfg.design==3)=2;
stat13 = ft_timelockstatistics(cfg, data13);

data42 = data;
x=find(conditions==4 | conditions==2);
data42.trial=data.trial(x,:,:);
conditions42 = conditions(x);
cfg.design = conditions42;
cfg.design(cfg.design==4)=1;
stat42 = ft_timelockstatistics(cfg, data42);

% hemiflip condition 4 and 2, and pool with conditions 1 and 3.
dataX = data;
x=find(conditions==4 | conditions==2);
tmp42 = Ci(x,:,:,:);
rev = [dim/2+1:dim 1:dim/2];
tmp42 = tmp42(:, rev, rev,:);

% Only take lower triangle of dim dimensions
tmp42 = reshape(tmp42, [size(tmp42,1), dim*dim, size(tmp42,4)]);
tmp42 = tmp42(:,sel,:);
tmp42 = reshape(tmp42, [size(tmp42,1), 1, size(tmp42,2)*size(tmp42,3)]);

dataX.trial(x,1,:) = tmp42;
cfg.design = conditions;
cfg.design(cfg.design==4)=1;
cfg.design(cfg.design==3)=2;
statX = ft_timelockstatistics(cfg, dataX);
%}

