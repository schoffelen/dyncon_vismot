function [statl,statr,sl,sr] = doSourcestatisticsCong2PreviousSmooth(band,correctm,nrand)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
load(['grandavgPrePreviousCong3Smooth-',band],'sl', 'sr');

if nargin<2, correctm = 'no'; end
if nargin<3, nrand    = 5000; end

nsubj = size(sl.pow,1);
npos  = size(sl.pos,1);
ntim  = size(sl.pow,4);

% global demeaning
%sl.pow = sl.pow ./ repmat(mean(sl.pow,4),[1 1 1 4]);
%sr.pow = sr.pow ./ repmat(mean(sr.pow,4),[1 1 1 4]);
%sl.pow = sl.pow - repmat(mean(sl.pow,1),[17 1 1 1]);
%sr.pow = sr.pow - repmat(mean(sr.pow,1),[17 1 1 1]);

%sl.pow(:,:,:,[1 3]) = sl.pow(:,:,:,[1 3])-repmat(mean(sl.pow(:,:,:,[1 3]),4),[1 1 1 2]);
%sl.pow(:,:,:,[2 4]) = sl.pow(:,:,:,[2 4])-repmat(mean(sl.pow(:,:,:,[2 4]),4),[1 1 1 2]);
%sr.pow(:,:,:,[1 3]) = sr.pow(:,:,:,[1 3])-repmat(mean(sr.pow(:,:,:,[1 3]),4),[1 1 1 2]);
%sr.pow(:,:,:,[2 4]) = sr.pow(:,:,:,[2 4])-repmat(mean(sr.pow(:,:,:,[2 4]),4),[1 1 1 2]);

sl.pow  = reshape(permute(mean(sl.pow,3), [1 4 2 3]), [nsubj*ntim npos]);
sl.powdimord = 'rpt_pos';
sl.freq = mean(sl.freq);
sr.pow  = reshape(permute(mean(sr.pow,3), [1 4 2 3]), [nsubj*ntim npos]);
sr.powdimord = 'rpt_pos';
sr.freq = mean(sr.freq);

sl.pow = double(sl.pow);
sr.pow = double(sr.pow);
sl = rmfield(sl, 'time');
sr = rmfield(sr, 'time');

design = [repmat(1:nsubj, [1 4]); ones(1,nsubj) ones(1,nsubj)*2 ones(1,nsubj) ones(1,nsubj)*2;ones(1,2*nsubj) ones(1,2*nsubj)*2];

cfg = [];
cfg.parameter   = 'pow';
cfg.design      = design;
cfg.ivar        = [2 3];
cfg.uvar        = 1;
cfg.statistic   = 'anova2x2rm';
cfg.method      = 'montecarlo';
cfg.numrandomization = nrand;
cfg.correctm    = correctm;
cfg.implementation = 'new';
cfg.resampling   = 'bootstrap';
cfg.precondition = 'after';
cfg.alpha = 0.025;
%statc = ft_sourcestatistics(cfg, sl);
%stati = ft_sourcestatistics(cfg, sr);
cfg.clusterthreshold = 'nonparametric_individual';
cfg.clusteralpha = 0.025;
dim = sl.dim;
for k = 1:size(sl.pow,1)
  fprintf('smoothing data\n');
  tmp = reshape(sl.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  %tmp(sl.inside) = tmp(sl.inside)-nanmedian(tmp(sl.inside));
  sl.pow(k,:) = tmp(:)';
  tmp2 = reshape(sr.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  %tmp2(sr.inside) = tmp2(sr.inside)-nanmedian(tmp2(sr.inside));
  sr.pow(k,:) = tmp2(:)';
end
cfg.f = 'main1';
statl(1) = ft_sourcestatistics(cfg, sl);
statr(1) = ft_sourcestatistics(cfg, sr);
cfg.f = 'main2';
statl(2) = ft_sourcestatistics(cfg, sl);
statr(2) = ft_sourcestatistics(cfg, sr);
cfg.f = 'interaction';
statl(3) = ft_sourcestatistics(cfg, sl);
statr(3) = ft_sourcestatistics(cfg, sr);

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
save(['statPrePreviousCong3Smooth-',band],'statl','statr');
