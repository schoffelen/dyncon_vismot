load list
for subj=18:19
  subj
  subjectname = list{subj};
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
  
  filename = sprintf('/project/3011085.03/analysis/artifact/%s_sourcepow.mat', subjectname);
  if isfile(filename)
    % add sampleinfo to data
    tmp = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
    data.sampleinfo = tmp.sampleinfo; clear tmp
    load(filename)
    tmpidx=[];
    % ft_rejectartifact doesn't work on this data, remove trials manually
    for l=1:size(artifact.summary.artifact,1)
      try
        tmpidx(l) = find(data.sampleinfo(:,1)==artifact.summary.artifact(l,1));
      end
    end
    data.trial(tmpidx,:,:) = [];
    data.sampleinfo(tmpidx,:)=[];
    data.trialinfo(tmpidx, :)=[];
  else
    cfg=[];
    cfg.method = 'summary';
    data = ft_rejectvisual(cfg, data);
    try
      artifact = data.cfg.artfctdef;
    catch
      artifact = tmpdata.cfg.previous.artfctdef;
    end
    save(filename, 'artifact')
  end
  
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
  
  idx1 = find(conditions==1);
  idx2 = find(conditions==2);
  idx3 = find(conditions==3);
  idx4 = find(conditions==4);
  
  % first input into SVM will be -1, second +1
  data31 = data;
  data31.trial = data.trial([idx3;idx1],:,:);
  data24 = data;
  data24.trial = data.trial([idx2;idx4],:,:);
  
  n=size(data.trial,1)/4;
  
  cfg=[];
  cfg.method= 'crossvalidate';
  cfg.mva=  {dml.standardizer dml.svm};
  cfg.statistic = {'confusion'  'accuracy'};
  cfg.type= 'nfold';
  cfg.nfolds = 5;
  cfg.design = [ones(n,1); 2*ones(n,1)];
  stat13 = ft_timelockstatistics(cfg, data31);
  stat13.statistic
  
  stat42 = ft_timelockstatistics(cfg, data24);
  stat42.statistic
   
  filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_svm1342.mat']);
  
  save(filename, 'stat13','stat42', 'cfg', 'data31', 'data24')
  
  % randomize
  numrandomization = 100;
  for k=1:numrandomization
    randseq = randperm(numel(cfg.design));
    
    qsubfeval(@vismot_execute_pipeline, 'vismot_bf_svm_rand', subject.name, {'randnr', k}, {'randseq', randseq}, 'memreq', 10*1024^3, 'timreq', 1800, 'batchid', sprintf('pow_svm_%s_rand%d', subjectname, k));
  end
  cd tmp
  keep subj list
end