function [source] = computeVoxeldataTFRavgZscore(subject, frequency);

[source] = computeFiltersTFR(subject, frequency);

cd('/analyse/1/Project0002/tmpProject0030');
load([subject.name,'tfrperiHanning2fourier']);

ninside = numel(inside);
npos    = size(source.pos,3);
ntim    = numel(allfreq{1}.time);

filt = source.filter;
if size(filt{inside(1)},1)>1,
  error('only vectorized filters are allowed');
end

filt = cat(1,filt{inside});

ix = nearest(allfreq{1}.freq, frequency);
for k = 1:4
  fourier = allfreq{k}.fourierspctrm(:,:,ix,:);
  siz     = size(fourier);
  fourier = reshape(fourier, [siz(1:2) siz(4)]);
  fourier = fourier./max(abs(fourier(:)));

  sumval  = zeros(ninside, ntim);
  ssqval  = zeros(ninside, ntim);
  for m = 1:ntim
    m
    tmp = transpose(fourier(:,:,m));
    tmp(~isfinite(tmp)) = 0;
    tmp2 = abs(filt*tmp).^2;
    sumval(:,m) = sum(tmp2,   2);
    ssqval(:,m) = sum(tmp2.^2,2);
  end
end

n = reshape(sum(isfinite(fourier(:,1,:)),1),[1 ntim]);
c = sumval./n(ones(1,ninside),:);
v = (ssqval - (sumval.^2)./n(ones(1,ninside),:))./(n(ones(1,ninside),:)-1);
