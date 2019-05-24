load list
n=numel(list);

numrandomization = 100;
nfreq = 4;
nfolds = 5;
load /project/3011085.03/scripts/fieldtrip/template/sourcemodel/standard_sourcemodel3d4mm.mat
ninsidepos = sum(sourcemodel.inside);

% prepare matrices
randprimal = zeros(numrandomization,6,ninsidepos,nfreq);
primal = zeros(n,6,ninsidepos,nfreq);

randacc = zeros(n,numrandomization);
for k=1:n
    k
    subjectname = list{k};
    subject = vismot_subjinfo(subjectname);
    filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_svm.mat']);
    M = load(filename, 'stat');
    
    for l=1:nfolds
        for m=1:6
            tmpprimal(l,m,:) = M.stat.model{l}.model{m}.primal;
        end
    end
    primal(k,:,:,:) = reshape(mean(tmpprimal,1), [6, ninsidepos, nfreq]);
    
    tmp2 = zeros(numrandomization,6, ninsidepos,nfreq);
    for l=1:numrandomization
        %       filename = fullfile(subject.pathname,'pow','rand', [subject.name, sprintf('_source3d4mm_pre_naive_rand%d.mat', l)]);
        filename = sprintf('/project_ext/3010029/reproducescript/vismsot/analysis/pow/rand/%s_source3d4mm_pre_svm_rand%d.mat', subject.name, l);
        randtmp = load(filename);
        
        randacc(k,l) = randtmp.stat.statistic.accuracy;
        tmp1 = zeros(nfolds, 6, ninsidepos, nfreq);
        
        for m=1:nfolds % number of folds
            for jj=1:6
            tmp1(m,jj,:,:) = reshape(randtmp.stat.model{m}.model{jj}.primal, ninsidepos, nfreq);
            end
        end
        tmp2(l,:,:,:) = squeeze(nanmean(tmp1,1));
    end
    randprimal(k, :, :, :) = nanmean(tmp2,1);

    acc(k,1) = M.stat.statistic.accuracy;
%     model{k} = M.model;
    p(k) = (sum(randacc(k,:)>=acc(k))+1)/(numrandomization+1); % See Ojala & Garriga 2010 (JMLR)
    
    diff{k} = sourcemodel;
    diff{k}.primal = zeros(6, numel(sourcemodel.inside), 4);
    diff{k}.primal(:,diff{k}.inside,:)=squeeze(primal(k,:,:,:))./squeeze(mean(randprimal,1));
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
    diff2{k}.primal = zeros(6,numel(sourcemodel.inside), 4);
    diff2{k}.primal(:,sourcemodel.inside,:) = squeeze(primal(k,:,:,:));
    diff2{k}.pos = sourcemodel.pos;
    nul{k} = diff2{k};
    nul{k}.primal = 0*nul{k}.primal;
end

cfg=[];
cfg.method = 'analytic';

cfg.statistic = 'depsamplesT';

cfg.correctm = 'no';
cfg.parameter = 'primal';
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






