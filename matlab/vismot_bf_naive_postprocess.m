load list
n=numel(list);

numrandomization = 500;
nfreq = 4;
nfolds = 5;
load /project/3011085.03/scripts/fieldtrip/template/sourcemodel/standard_sourcemodel3d4mm.mat
ninsidepos = sum(sourcemodel.inside);

% prepare matrices
randdivergence = zeros(n,ninsidepos,nfreq);
randacc = zeros(n,numrandomization);
for k=1:n
    k
    subjectname = list{k};
    subject = vismot_subjinfo(subjectname);
    filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_naive.mat']);
    M = load(filename, 'stat','model', 'stat_perfreq');
    
    tmp2 = zeros(numrandomization,ninsidepos,nfreq);
    for l=1:numrandomization
        if l==100 || l==200 || l==300 || l==400
            l
        end
        %       filename = fullfile(subject.pathname,'pow','rand', [subject.name, sprintf('_source3d4mm_pre_naive_rand%d.mat', l)]);
        filename = sprintf('/project_ext/3010029/reproducescript/vismsot/analysis/pow/rand/%s_source3d4mm_pre_naive_rand%d.mat', subject.name, l);
        randtmp = load(filename);
        
        randacc(k,l) = randtmp.stat.statistic.accuracy;
        tmp1 = zeros(nfolds, ninsidepos, nfreq);
        
        for m=1:nfolds % number of folds
            tmp1(m,:,:) = reshape(randtmp.stat.model{m}.divergence, ninsidepos, nfreq);
        end
        tmp2(l,:,:) = nanmean(tmp1,1);
    end
    randdivergence(k, :, :) = nanmean(tmp2,1);

    acc(k,1) = M.stat.statistic.accuracy;
    model{k} = M.model;
    p(k) = (sum(randacc(k,:)>=acc(k))+1)/(numrandomization+1); % See Ojala & Garriga 2010 (JMLR)
    
    diff{k} = model{k};
    diff{k}.divergence(diff{k}.inside,:)=model{k}.divergence(diff{k}.inside,:)./squeeze(randdivergence(k,:,:));
end


figure;
for k=1:n
    subplot(4,5,k);hold on
    histogram(randacc(k,:)); vline(acc(k));
    title(sprintf('p=%d', p(k)));
end


idx_sign = find(p<0.05);
nsign = numel(idx_sign);
diff2 = diff(idx_sign);

for k=1:nsign
    diff2{k}.pos = diff2{1}.pos;
    nul{k} = diff2{k};
    nul{k}.divergence = 0*nul{k}.divergence;
end

cfg=[];
cfg.method = 'analytic';

cfg.statistic = 'depsamplesT';

cfg.correctm = 'no';
cfg.parameter = 'divergence';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,nsign) 2*ones(1,nsign); 1:nsign, 1:nsign];

% cfg.method = 'montecarlo';
% cfg.numrandomization = 1000;
% cfg.correctm = 'cluster';
stat = ft_sourcestatistics(cfg, diff2{:}, nul{:})

cfg2 = [];
% cfg2.frequency = 10;
cfg2.parameter = {'stat', 'mask'};
stat_int = ft_sourceinterpolate(cfg2, stat, mri)






