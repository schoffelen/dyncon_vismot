function [freq] = checkFreqpre(freqpre, savename)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');

if nargin<2,
  savename = '';
end

if isstr(freqpre),
  load(freqpre);
end

if iscell(freqpre), tmp = selectdata(freqpre{:}, 'param', 'fourierspctrm'); end
for k = 1:numel(tmp.freq)
  k
  tmpx = selectdata(tmp, 'foilim', tmp.freq(k).*[1 1]);
  tmpx = checkdata(tmpx, 'cmbrepresentation', 'full');
  tmpx.crsspctrm = tmpx.crsspctrm./max(tmpx.crsspctrm(:));
  for m = 1:size(tmpx.crsspctrm,1)
    nc(m,k) = norm(squeeze(tmpx.crsspctrm(m,:,:)),'fro');
  end
end
nc   = standardise(nc,1);
snc  = sum(nc>2,2);
sel  = find(snc<3);
freq = selectdata(tmp, 'rpt', sel);

if ~isempty(savename)
  save(savename, 'freq');
end
