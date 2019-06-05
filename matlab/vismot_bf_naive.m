addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

subject = vismot_subjinfo(subjectname);
freqs = [10 20 40 60];
k=1;
for f=freqs
  filename = fullfile(subject.pathname,'pow', [subject.name,sprintf('_source3d4mm_%s_', 'pre'), num2str(f,'%03d')]);
  S{k} = load(filename, 'sourcepow');
  k=k+1;
end

ntrials = size(S{1}.sourcepow.pow,1);
ninsidepos = sum(S{1}.sourcepow.inside);
nfreq = numel(freqs);
insidepow = zeros(ntrials, ninsidepos, nfreq);
for k=1:numel(S)
  insidepow(:,:,k) = S{k}.sourcepow.pow(:,S{k}.sourcepow.inside);
end
data.trial = reshape(insidepow, [ntrials, 1, ninsidepos*nfreq]);
data.time = 1:size(data.trial,3);
data.dimord = 'rpt_chan_time';
data.trialinfo = S{1}.sourcepow.trialinfo;
data.label = {'pow'};

idx = find(data.trialinfo(:,end)==5);
data.trial(idx,:,:)=[];
data.trialinfo(idx,:)=[];
insidepow(idx,:,:) = [];

cfg=[];
cfg.method = 'summary';
data = ft_rejectvisual(cfg, data);

% equalize conditions
ntrials = min([sum(data.trialinfo(:,end)==1), sum(data.trialinfo(:,end)==2), sum(data.trialinfo(:,end)==3), sum(data.trialinfo(:,end)==4)]);
tmpidx = [];
for k=1:numel(unique(data.trialinfo(:,end)))
  tmp = find(data.trialinfo(:,end)==k);
  P = randperm(numel(tmp));
  tmpidx = [tmpidx; tmp(P(1:ntrials))];
end
data.trial = data.trial(tmpidx, :, :);
data.trialinfo = data.trialinfo(tmpidx,:);

conditions = data.trialinfo(:, end);
insidepow = insidepow(tmpidx, :,:);

cfg=[];
cfg.method= 'crossvalidate';
cfg.mva=  {dml.standardizer dml.naive};%dml.one_against_one('mva', dml.svm)}
cfg.statistic = {'confusion'  'accuracy'};
cfg.type= 'nfold';
cfg.nfolds = 5;
cfg.design = conditions;
stat = ft_timelockstatistics(cfg, data);
stat.statistic

model = removefields(S{1}.sourcepow, {'pow', 'dimord', 'cfg', 'trialinfo'});
for k=1:cfg.nfolds
divergence(k,:,:) = reshape(stat.model{k}.divergence, ninsidepos, 4);
end
model.freq = freqs;
model.divergence = zeros(numel(model.inside),nfreq);
model.divergence(model.inside,:) = squeeze(nanmean(divergence,1));



%% seperately for each frequency
for k=1:nfreq
  k
  data_perfreq{k} = removefields(data, 'trial');
  data_perfreq{k}.trial = permute(insidepow(:,:,k), [1, 3, 2]);
  data_perfreq{k}.time = 1:ninsidepos;
  stat_perfreq{k} = ft_timelockstatistics(cfg, data_perfreq{k});
end

filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_naive.mat']);
% filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_svm.mat']);

save(filename, 'stat','model', 'stat_perfreq', 'cfg', 'data')





% randomize
numrandomization = 500;
for k=1:numrandomization
  randseq = randperm(numel(cfg.design));
%   vismot_execute_pipeline('vismot_bf_naive_rand', subject.name, {'randnr', k}, {'randseq', randseq},{'data', data}, {'cfg',cfg});

  qsubfeval(@vismot_execute_pipeline, 'vismot_bf_naive_rand', subject.name, {'randnr', k}, {'randseq', randseq}, 'memreq', 10*1024^3, 'timreq', 1800, 'batchid', sprintf('pow_naive_%s_rand%d', subjectname, k));
end


%{

randacc = zeros(numrandomization,1);
for k=1:numrandomization
  k
  cfg.design = cfg.design(randperm(numel(cfg.design)));
  tmpstat{k} = ft_timelockstatistics(cfg, data);
  randacc(k,1) = tmpstat{k}.statistic.accuracy;
end
%}