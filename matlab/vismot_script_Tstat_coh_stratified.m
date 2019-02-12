datadir = '/project/3011085.03/analysis/source/lambda/';
foi=10;
lambda='10%';
N=19;

for k=1:N
  tmp{k} = load(fullfile([datadir sprintf('%s/%s_coh6d10mm_%03d.mat', lambda, list{k}, foi)]));
end
sourcemodel = tmp{1}.sourcemodel;

source = sourcemodel;
source.coh13 = nan(N, size(sourcemodel.inside,1), size(sourcemodel.inside,1));
source.coh42 = nan(N, size(sourcemodel.inside,1), size(sourcemodel.inside,1));
insidepos = find(sourcemodel.inside);
for k=1:19
  source.coh13(k,insidepos,insidepos) = tmp{k}.zx13;
  source.coh42(k,insidepos,insidepos) = tmp{k}.zx42;
end

source.frequency = foi;
source.dimord = 'rpt_pos_time';
source.time = 1:6804;

nul = source;
nul.coh13(:) = 0;
nul.coh42(:) = 0;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,N) ones(1,N)*2;1:N 1:N];
cfgs.correctm = 'no';
cfgs.numrandomization = 1000;
cfgs.clusteralpha = 0.05;
cfgs.correcttail = 'prob';
cfgs.parameter =  'coh13';
stat13 = ft_sourcestatistics(cfgs, source, nul);