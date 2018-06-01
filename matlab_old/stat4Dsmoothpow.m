subjinfo;
names = {SUBJ(:).name};
names = names([1:3 5 6 8:11 13:16 18:20]);

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);
for k = 1:numel(names)
  k
  fname = [names{k},'powFastSlowSmooth'];
  load(fname);
  if k == 1,
    dat = zeros(prod(s{1}.dim),numel(names));
    dim = s{1}.dim;
  end
  dum = reshape(mean(s{1}.stat+s{2}.stat+s{3}.stat+s{4}.stat,2)./4,dim);
  spm_smooth(dum,dum,1.5);
  dat(:,k) = dum(:);
  clear dum;
  clear s;
end
load(fname);
sd        = s{1};clear s;
sd.pow    = dat;
sd.pow(:,17:32) = ones(size(sd.pos,1),1)*nanmean(sd.pow(sd.inside,1:16),1);
sd.powdimord = 'pos_rpt_freq';
sd.freq   = 10;
sd = rmfield(sd, 'stat');
sd = rmfield(sd, 'prob');
sd = rmfield(sd, 'mask');
sd = rmfield(sd, 'fwhm');

cfg = [];
cfg.implementation = 'new';
cfg.method = 'montecarlo';
cfg.numrandomization = 1000;
cfg.statistic      = 'pooledT';
cfg.parameter      = 'pow';
cfg.correctm       = 'no';
cfg.design         = [ones(1,16) ones(1,16)*2; 1:16 1:16];
cfg.ivar           = 1;
cfg.uvar           = 2;
stat = ft_sourcestatistics(cfg, sd);

