function [statc,stati,statx] = doSourcestatisticsCong2PreviousSmooth(band,correctm,nrand)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
%load(['grandavgPrePreviousCong4Smooth-',band],'sc', 'si');
load(['grandavgPrePreviousCong5',band],'sc', 'si');

if nargin<2, correctm = 'no'; end
if nargin<3, nrand = 5000; end

nsubj = size(sc.pow,1);

sc.pow = mean(sc.pow,3);
sc.freq = mean(sc.freq);
sc.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(sc.pow,2),[1 size(sc.pow,2) 1]);
si.pow = mean(si.pow,3);
si.freq = mean(si.freq);
si.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(si.pow,2),[1 size(si.pow,2) 1]);


cfg = [];
cfg.parameter   = 'pow';
cfg.design      = [ones(1,nsubj) 2*ones(1,nsubj);1:nsubj 1:nsubj];
cfg.ivar        = 1;
cfg.uvar        = 2;
%cfg.statistic   = 'pooledTtfce';
%cfg.statistic   = 'pooledT';
cfg.statistic   = 'depsamplesT';
%cfg.dim    = si.dim;
%cfg.inside = si.inside;
cfg.method      = 'montecarlo';
cfg.numrandomization = nrand;
cfg.correctm    = correctm;
cfg.implementation = 'new';
cfg.alpha = 0.025;
cfg.dim = sc.dim;
%statc = ft_sourcestatistics(cfg, sc);
%stati = ft_sourcestatistics(cfg, si);
cfg.clusterthreshold = 'nonparametric_individual';
%cfg.clusterstatistic = 'maxsize';
cfg.clusteralpha = 0.05;

dim = sc.dim;
sx  = sc;
for k = 1:nsubj
  tmp = reshape(sc.pow(k,:),dim);
  tmp2 = reshape(si.pow(k,:),dim);
  
  tmp3 = (tmp+flipdim(tmp2,1))./sqrt(2);

  spm_smooth(tmp,tmp,2);
  spm_smooth(tmp2,tmp2,2);
  spm_smooth(tmp3,tmp3,2);
  sc.pow(k,:) = tmp(:)';
  sc.pow(k+nsubj,:) = nanmean(sc.pow(k,:));
  %sc.pow(k+nsubj,:) = 0;
  si.pow(k,:) = tmp2(:)';
  si.pow(k+nsubj,:) = nanmean(si.pow(k,:));
  %si.pow(k+nsubj,:) = 0;
  sx.pow(k,:) = tmp3(:)'; 
  %sx.pow(k+nsubj,:) = 0;
  sx.pow(k+nsubj,:) = nanmean(sx.pow(k,:));
end
cfg.numrandomization = 1
cfg.correctm = 'no';
statc = ft_sourcestatistics(cfg, sc);
stati = ft_sourcestatistics(cfg, si);
cfg.numrandomization = nrand;
cfg.correctm = correctm;
statx = ft_sourcestatistics(cfg, sx);

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
%save(['statPrePreviousCong4Smooth-',band],'statc','stati','statx');
save(['statPrePreviousCong5-',band],'statc','stati','statx');
