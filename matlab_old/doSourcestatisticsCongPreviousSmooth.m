function [statlc,statli,statrc,statri,statv1,statv2,statm1,statm2,statv] = doSourcestatisticsCongPreviousSmooth(band,correctm)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/powprepreviouscong20101210');
load(['grandavgPrePreviousCongSmooth-',band],'slc', 'sli','src','sri');

if nargin<2, correctm = 'no'; end

nsubj = size(slc.pow,1);

slc.pow = mean(slc.pow,3);
slc.freq = mean(slc.freq);
slc.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(slc.pow,2),[1 size(slc.pow,2) 1]);
sli.pow = mean(sli.pow,3);
sli.freq = mean(sli.freq);
sli.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(sli.pow,2),[1 size(sli.pow,2) 1]);
src.pow = mean(src.pow,3);
src.freq = mean(src.freq);
src.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(src.pow,2),[1 size(src.pow,2) 1]);
sri.pow = mean(sri.pow,3);
sri.freq = mean(sri.freq);
sri.pow((nsubj+1):(nsubj*2),:,:) = repmat(nanmean(sri.pow,2),[1 size(sri.pow,2) 1]);


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
statlc = ft_sourcestatistics(cfg, slc);
statli = ft_sourcestatistics(cfg, sli);
statrc = ft_sourcestatistics(cfg, src);
statri = ft_sourcestatistics(cfg, sri);

% what could be interesting is based on the following argumentation
% if the next trial is (in)congruent, there is a behavioral benefit if the previous
% trial was (in)congruent.
% Thus, for regions processing the target stimulus lc - ri, and rc - li are interesting
% Thus, for regions more related to motor execution lc - li, and rc - ri are interesting
% in order to capture the behaviorally relevant regions

dim = statlc.dim;
s   = slc;
for k = 1:nsubj
  tmp = reshape(slc.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  tmp2 = reshape(sri.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  s.pow(k,:) = (tmp(:)-tmp2(:))./sqrt(2);
  s.pow(k+nsubj,:) = nanmean(s.pow(k,:));
end
inside = zeros(dim);
inside(s.inside) = 1;
statv1 = ft_sourcestatistics(cfg, s);

dim = statlc.dim;
s2   = src;
for k = 1:nsubj
  tmp = reshape(src.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  tmp2 = reshape(sli.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  s2.pow(k,:) = (tmp(:)-tmp2(:))./sqrt(2);
  s2.pow(k+nsubj,:) = nanmean(s2.pow(k,:));
end
inside = zeros(dim);
statv2 = ft_sourcestatistics(cfg, s2);

for k = 1:nsubj*2
  tmp = reshape(s2.pow(k,:),dim);
  tmp = flipdim(tmp, 1);
  s.pow(k,:) = (s.pow(k,:) + tmp(:)')./sqrt(2);
end
statv = ft_sourcestatistics(cfg, s);

dim = statlc.dim;
s   = slc;
for k = 1:nsubj
  tmp = reshape(slc.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  tmp2 = reshape(sli.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  s.pow(k,:) = (tmp(:)-tmp2(:))./sqrt(2);
  s.pow(k+nsubj,:) = nanmean(s.pow(k,:));
end
inside = zeros(dim);
inside(s.inside) = 1;
statm1 = ft_sourcestatistics(cfg, s);

dim = statlc.dim;
s   = src;
for k = 1:nsubj
  tmp = reshape(src.pow(k,:),dim);
  spm_smooth(tmp,tmp,1.5);
  tmp2 = reshape(sri.pow(k,:),dim);
  spm_smooth(tmp2,tmp2,1.5);
  s.pow(k,:) = (tmp(:)-tmp2(:))./sqrt(2);
  s.pow(k+nsubj,:) = nanmean(s.pow(k,:));
end
inside = zeros(dim);
statm2 = ft_sourcestatistics(cfg, s);

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/powprepreviouscong20101210');
save(['statPrePreviousCongSmooth-',band],'statlc','statli','statrc','statri','statv1','statv2','statm1','statm2','statv');
