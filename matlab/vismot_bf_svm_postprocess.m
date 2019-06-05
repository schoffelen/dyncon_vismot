load list
n=numel(list);

numrandomization = 100;
nfreq = 4;
nfolds = 5;
load /project/3011085.03/scripts/fieldtrip/template/sourcemodel/standard_sourcemodel3d4mm.mat
ninsidepos = sum(sourcemodel.inside);
freqs = [10 20 40 60];
% prepare matrices
% primal = zeros(n,6,ninsidepos,nfreq);
primal = zeros(n,ninsidepos,nfreq);

randacc = zeros(n,numrandomization);
for q=[1 2];
for k=1:n
    k
    subjectname = list{k};
    subject = vismot_subjinfo(subjectname);
    filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_svm1342.mat']);
    M = load(filename);
    if q==1
      M.stat=M.stat13;
    elseif q==2
      M.stat = M.stat42;
    end
    for l=1:nfolds
%         for m=1%:6
%             tmpprimal(l,m,:) = M.stat.model{l}.model{m}.primalinv;
            tmpprimal(l,:) = M.stat.model{l}.primal;
%         end
    end
%     primal(k,:,:,:) = reshape(mean(tmpprimal,1), [6, ninsidepos, nfreq]);
if q==1
        primal13(k,:,:) = reshape(mean(tmpprimal,1), [ninsidepos, nfreq]);
elseif q==2
  primal42(k,:,:) = reshape(mean(tmpprimal,1), [ninsidepos, nfreq]);
end

    for l=1:numrandomization
        %       filename = fullfile(subject.pathname,'pow','rand', [subject.name, sprintf('_source3d4mm_pre_naive_rand%d.mat', l)]);
        filename = sprintf('/project_ext/3010029/reproducescript/vismsot/analysis/pow/rand/%s_source3d4mm_pre_svm1342_rand%d.mat', subject.name, l);
        randtmp = load(filename);
        if q==1
          randtmp.statistic = randtmp.statistic13;
          randacc13(k,l) = randtmp.statistic.accuracy;
        elseif q==2
          randtmp.statistic = randtmp.statistic42;
          randacc42(k,l) = randtmp.statistic.accuracy;
        end
    end
if q==1
    acc13(k,1) = M.stat.statistic.accuracy;
elseif q==2;
  acc42(k,1) = M.stat.statistic.accuracy;
end
%     model{k} = M.model;
if q==1
      p(k) = (sum(randacc13(k,:)>=acc13(k))+1)/(numrandomization+1); % See Ojala & Garriga 2010 (JMLR)
elseif q==2
  p(k) = (sum(randacc42(k,:)>=acc42(k))+1)/(numrandomization+1); % See Ojala & Garriga 2010 (JMLR)
end
end
end

figure;
for k=1:n
    subplot(4,5,k);hold on
    histogram(randacc42(k,:)); vline(acc42(k));
    title(sprintf('p=%d', p(k)));
end

for k=1:n
weights{k} = sourcemodel;
weights{k}.stat = zeros(numel(sourcemodel.inside), nfreq);
% weights{k}.stat(sourcemodel.inside,:) = squeeze(squeeze(nanmean(primal(k,5,:,:),2)));
weights{k}.stat(sourcemodel.inside,:) = squeeze(primal(k,:,:));
weights{k}.dimord = 'pos_freq';
weights{k}.freq = freqs;

nul{k} = weights{k};
nul{k}.stat = 0*nul{k}.stat;
end

% W = weights{1};
% W.stat(W.inside,:) = squeeze(squeeze(nanmean(nanmean(primal,2),1)));


cfg=[];
cfg.method = 'analytic';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'no';
cfg.parameter = 'stat';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,n) 2*ones(1,n); 1:n, 1:n];
stat = ft_sourcestatistics(cfg, weights{:}, nul{:});

cfgp=[];
cfgp.funparameter = 'stat';
cmap = flipud(brewermap(64,'RdBu'));
cfgp.funcolormap = cmap;
cfgp.maskparameter = 'stat';
cfgp.method = 'slice';
for k=freqs
  cfgp.frequency=k;
  ft_sourceplot(cfgp, stat);
end
