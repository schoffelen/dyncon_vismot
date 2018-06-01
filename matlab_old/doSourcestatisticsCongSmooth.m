function [stat,stat13,stat42] = doSourcestatisticsCongSmooth(band,correctm)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source');
load(['grandavgCongSmooth-',band],'s13', 's42');

if nargin<2, correctm = 'no'; end

nsubj = size(s13.pow,1);

s13.pow = mean(s13.pow,3);
s13.freq = mean(s13.freq);
s13.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(s13.pow,2),[1 size(s13.pow,2) 1]);
s42.pow = mean(s42.pow,3);
s42.freq = mean(s42.freq);
s42.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(s42.pow,2),[1 size(s42.pow,2) 1]);


cfg = [];
cfg.parameter   = 'pow';
cfg.design      = [ones(1,nsubj) 2*ones(1,nsubj);1:nsubj 1:nsubj];
cfg.ivar        = 1;
cfg.uvar        = 2;
cfg.statistic   = 'pooledT';
cfg.method      = 'montecarlo';
cfg.numrandomization = 5000;
cfg.correctm    = correctm;
cfg.implementation = 'new';
cfg.alpha = 0.025;
stat13 = ft_sourcestatistics(cfg, s13);
stat42 = ft_sourcestatistics(cfg, s42);

% do flipping and smoothing
dim = stat13.dim;
s   = s13;
for k = 1:nsubj*2
  tmp = reshape(s13.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  tmp  = flipdim(tmp, 1);
  tmp2 = reshape(s42.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  s.pow(k,:) = (tmp(:)+tmp2(:))./sqrt(2);
end
inside = zeros(dim);
inside(s.inside) = 1;
inside = flipdim(inside,1);
inside(s.inside) = inside(s.inside) + 1;
s.inside = find(inside==2);
stat = ft_sourcestatistics(cfg, s);

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source');
save(['statCongSmooth-',band],'stat13', 'stat42','stat');
