function [stat13, stat42] = doFreqstatisticsTFR(subject, band)

cd(subject.pathname);
cd('freq/tfr');
load([subject.name,'tfr',band]);

tim1 = max(freq{1}.time(1),freq{3}.time(1));
tim2 = min(freq{1}.time(end),freq{3}.time(end));
tim3 = max(freq{2}.time(1),freq{4}.time(1));
tim4 = min(freq{2}.time(end),freq{4}.time(end));
tim1 = max(tim1,tim3);
tim2 = min(tim2,tim4);

freq{1} = checkdata(freq{1}, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
freq{2} = checkdata(freq{2}, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
freq{3} = checkdata(freq{3}, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
freq{4} = checkdata(freq{4}, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});

freq{1} = selectdata(freq{1}, 'toilim', [tim1 tim2]);
freq{2} = selectdata(freq{2}, 'toilim', [tim1 tim2]);
freq{3} = selectdata(freq{3}, 'toilim', [tim1 tim2]);
freq{4} = selectdata(freq{4}, 'toilim', [tim1 tim2]);

tmp1           = selectdata(freq{1},freq{3},'param','powspctrm');
tmp1.rt        = [freq{1}.rt(:);freq{3}.rt(:)];
tmp1.condition = [ones(numel(freq{1}.cumtapcnt),1); -ones(numel(freq{3}.cumtapcnt),1)];
tmp1.cfg.trl   = zeros(numel(tmp1.rt),3);
tmp1.powspctrm = log10(tmp1.powspctrm);

tmp2           = selectdata(freq{4},freq{2},'param','powspctrm');
tmp2.rt        = [freq{4}.rt(:);freq{2}.rt(:)];
tmp2.condition = [ones(numel(freq{4}.cumtapcnt),1); -ones(numel(freq{2}.cumtapcnt),1)];
tmp2.cfg.trl   = zeros(numel(tmp2.rt),3);
tmp2.powspctrm = log10(tmp2.powspctrm);

cfg               = [];
cfg.method        = 'montecarlo';
cfg.statistic     = 'glm';
cfg.numrandomization = 0;
cfg.glm.statistic = 'T';
design1           = [tmp1.condition(:)';tmp1.rt(:)'];
design2           = [tmp2.condition(:)';tmp2.rt(:)'];
warning off;
for k = 1:numel(tmp1.time)
  tmpx      = selectdata(tmp1, 'toilim', tmp1.time(k));
  tmpdesign = design1;
  ix        = find(~isnan(tmpx.powspctrm(:,1)));
  tmpdesign = tmpdesign(:,ix);
  tmpx      = selectdata(tmpx, 'rpt', ix);
  condvec   = tmpdesign(1,:);
  n1        = sum(condvec == 1); ix1 = randperm(n1);
  n2        = sum(condvec ==-1); ix2 = randperm(n2);
  n         = min(n1,n2); ix1 = sort(ix1(1:n)); ix2 = sort(ix2(1:n));
  ix        = [ix1(:);ix2(:)+n1];
  tmpdesign = tmpdesign(:,ix);
  tmpx      = selectdata(tmpx, 'rpt', ix);
  cfg.design = orthogonalise(tmpdesign')';
  cfg.glm.contrast = [1 0];
  tmp13{k} = freqstatistics(cfg, tmpx);
  tmp13{k}.powspctrm = []; %hack to make selectdata work
  cfg.glm.contrast = [0 1];
  tmp13b{k} = freqstatistics(cfg, tmpx);
  tmp13b{k}.powspctrm = []; %hack to make selectdata work
  
  tmpx      = selectdata(tmp2, 'toilim', tmp1.time(k));
  tmpdesign = design2;
  ix        = find(~isnan(tmpx.powspctrm(:,1)));
  tmpdesign = tmpdesign(:,ix);
  condvec   = tmpdesign(1,:);
  tmpx      = selectdata(tmpx, 'rpt', ix);
  n1        = sum(condvec == 1); ix1 = randperm(n1);
  n2        = sum(condvec ==-1); ix2 = randperm(n2);
  n         = min(n1,n2); ix1 = sort(ix1(1:n)); ix2 = sort(ix2(1:n));
  ix        = [ix1(:);ix2(:)+n1];
  tmpdesign = tmpdesign(:,ix);
  tmpx      = selectdata(tmpx, 'rpt', ix);
  cfg.design = orthogonalise(tmpdesign')';
  cfg.glm.contrast = [1 0];
  tmp42{k} = freqstatistics(cfg, tmpx);
  tmp42{k}.powspctrm = []; %hack to make selectdata work
  cfg.glm.contrast = [0 1];
  tmp42b{k} = freqstatistics(cfg, tmpx);
  tmp42b{k}.powspctrm = []; %hack to make selectdata work
end 
stat13         = selectdata(tmp13{:},  'param', 'stat');
stat13.stat    = squeeze(mean(stat13.stat,2));
stat13.dimord  = 'chan_time';
stat13         = rmfield(stat13,  'freq');
stat13b        = selectdata(tmp13b{:}, 'param', 'stat');
stat13b.stat   = squeeze(mean(stat13b.stat,2));
stat13b.dimord = 'chan_time';
stat13b        = rmfield(stat13b, 'freq');
stat42         = selectdata(tmp42{:},  'param', 'stat');
stat42.stat    = squeeze(mean(stat42.stat,2));
stat42.dimord  = 'chan_time';
stat42         = rmfield(stat42,  'freq');
stat42b        = selectdata(tmp42b{:}, 'param', 'stat');
stat42b.stat   = squeeze(mean(stat42b.stat,2));
stat42b.dimord = 'chan_time';
stat42b        = rmfield(stat42b, 'freq');
keyboard

function [x] = orthogonalise(input)

input = standardise(input,1);
x     = input(:,1);

for i = 2:size(input,2)
  y   = input(:,i);
  y   = y-x*(pinv(x'*x)*x'*y);
  if any(y), x = [x y]; end
end
