cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
cd('powpstcong20101121/stats');
addpath('/home/jan/matlab/fieldtrip/external/spm2');
mri=ft_read_mri('/home/jan/matlab/mri/templateMRI.mnc');

for k = 1:4
  if k==1, 
    str = 'alpha';
  elseif k==2,
    str = 'beta';
  elseif k==3,
    str = 'gamma1';
  elseif k==4,
    str = 'gamma2';
  end
  load(['statCongSmooth-',str]);
  stat.prob(isnan(stat.prob)) = 1;
  stat13.prob(isnan(stat13.prob)) = 1;
  stat42.prob(isnan(stat42.prob)) = 1;
  stat.stat(isnan(stat.stat)) = 0;
  stat13.stat(isnan(stat13.stat)) = 0;
  stat42.stat(isnan(stat42.stat)) = 0;
  
  cfg           = [];
  cfg.parameter = {'prob' 'stat'};
  iall = ft_sourceinterpolate(cfg, stat,   mri);
  i13  = ft_sourceinterpolate(cfg, stat13, mri);
  i42  = ft_sourceinterpolate(cfg, stat42, mri);
  
  iall.stat(isnan(iall.stat)) = 0;
  i13.stat(isnan(i13.stat)) = 0;
  i42.stat(isnan(i42.stat)) = 0;
  iall.prob(isnan(iall.prob)) = 1;
  i13.prob(isnan(i13.prob)) = 1;
  i42.prob(isnan(i42.prob)) = 1;

  cfgp              = [];
  cfgp.method       = 'slice';
  cfgp.funparameter = 'stat';
  %cfgp.colorbar     = 'no';
  cfgp.nslices      = 16;
  cfgp.slicerange   = [25 155]
  cfgp.maskparameter = 'stat';
  figure;ft_sourceplot(cfgp, iall);title([str,'-pooled']);print(gcf, '-depsc2', [str,'-pooled']);
  figure;ft_sourceplot(cfgp, i13); title([str,'-left']);print(gcf, '-depsc2', [str,'-left']);
  figure;ft_sourceplot(cfgp, i42); title([str,'-right']);print(gcf, '-depsc2', [str,'-right']);
  cfgp.maskparameter = 'prob';
  cfgp.opacitylim    = [0.025 0.05];
  cfgp.opacitymap    = 'rampdown';
  figure;ft_sourceplot(cfgp, iall);title([str,'-pooledPmask']);print(gcf, '-depsc2', [str,'-pooledPmask']);
  figure;ft_sourceplot(cfgp, i13); title([str,'-leftPmask']);print(gcf, '-depsc2', [str,'-leftPmask']);
  figure;ft_sourceplot(cfgp, i42); title([str,'-rightPmask']);print(gcf, '-depsc2', [str,'-rightPmask']);

  clear iall i13 i42;
end
