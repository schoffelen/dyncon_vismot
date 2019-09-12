load vismot_parcels
cd /project/3011085.03/analysis/mim
load list

subjectname = list{1};

x = load(sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subjectname, 'previous'));
T = numel(unique(x.trialnumber));

filename = sprintf('%s_mim_pre_all', subjectname);
mim_all = load(filename);
c = mim_all.mim.mimspctrm;

filename = ['singletrial/', subjectname, '/', filename, '*'];
d = dir(filename);
dim = 16;

if dim~=6
  P3=1;
end
if dim~=16 && dim~=6
  P2=1;
end

freqs = 0.5:0.5:120;
nfreq = 159;%4;

% Compute pseudo values and average within regions/frequencies of interest.
Ci = zeros(numel(d), dim, dim, nfreq);
for k = 1:numel(d)
  k
  trialnumber(k) = str2num(char(extractBetween(d(k).name,'all_', '.mat')));
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
  elseif nfreq==80
    Ci(k,:,:,:) = Citmp2(:,:,3:2:2*(nfreq+1));
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
% save(filename, 'Ci', 'conditions', 'trialnumber', 'trialinfo');

% equalize conditions
ntrials = min([sum(conditions==1), sum(conditions==2), sum(conditions==3), sum(conditions==4), sum(conditions==5)]);
tmpidx = [];
for k=1:5
  tmp = find(conditions==k);
  P = randperm(sum(conditions==k));
  tmpidx = [tmpidx; tmp(P(1:ntrials))];
end
Ci = Ci(tmpidx,:,:,:);
conditions = conditions(tmpidx);

% Only take lower triangle of dim dimensions
if 0 % both if and else should lead to the same result. But not when reshape dimensions are wrong.
  Citmp = reshape(permute(Ci,[2,3,1,4]), [dim, dim, size(Ci,1)*nfreq]);
  Citmp2=zeros(size(Citmp));
  for k=1:size(Citmp,3)
    Citmp2(:,:,k) = tril(Citmp(:,:,k),-1);
  end
  
  % MAKE SURE DIMENSIONS ARE CORRECT. permute(reshape(Citmp2, [dim, dim,
  % nfreq, size(Ci,1)]), [4,1,2,3]); leads to accuracies of 50%!!
  Citmp = permute(reshape(Citmp2, [dim, dim, size(Ci,1), nfreq]), [3,1,2,4]); 
  Ci_small = single(Citmp);
  
  data.trial = reshape(Ci_small, [size(Ci,1), 1, dim*dim*nfreq]);
  
else
  sel = tril(ones(dim),-1)==1;
  Citmp = reshape(Ci, [size(Ci,1), dim*dim, size(Ci,4)]);
  Citmp = Citmp(:,sel,:);
  Citmp = reshape(Citmp, [size(Citmp,1), 1, size(Citmp,2)*size(Citmp,3)]);
  
  data.trial = Citmp;
end
data.time = 1:size(data.trial,3);
data.label = {'chan01'};
data.dimord = 'rpt_chan_time';
data.trialinfo = conditions(:);

% remove conditions 5
%{
x = find(conditions==5);
data.trial(x,:,:)=[];
conditions(x)=[];
%}

addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
cfg=[];
cfg.method = 'crossvalidate';
cfg.mva = {dml.naive};
cfg.statistic = {'confusion', 'accuracy'};
cfg.design = conditions;
cfg.type = 'nfold';
cfg.nfolds = 10;
cfg.resample = 0; % resamples conditions with fewer trials, and throws away oversampled conditions
stat = ft_timelockstatistics(cfg, data);
stat.statistic


numrandomization = 2;
r=zeros(numrandomization,1);
for k=1:numrandomization
  k
  cfg.design = conditions(randperm(numel(conditions)));
  randstat{k} = ft_timelockstatistics(cfg, ft_selectdata(cfgsel, data));
  r(k) = randstat{k}.statistic.accuracy;
end

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
tmp = reshape(permute(tmp42,[2,3,1,4]), [dim, dim, nfreq*size(tmp42,1)]);
for k=1:size(tmp,3)
  tmp(:,:,k) = tril(tmp(:,:,k));
end
tmp42 = permute(reshape(tmp, [dim, dim, nfreq, size(tmp42,1)]), [4,1,2,3]);
dataX.trial(x,1,:) = reshape(tmp42, [size(tmp42,1), 1, dim*dim*nfreq]);
cfg.design = conditions;
cfg.design(cfg.design==4)=1;
cfg.design(cfg.design==3)=2;
statX = ft_timelockstatistics(cfg, dataX);
%}

