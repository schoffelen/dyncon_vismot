datadir = '/project/3011085.03/analysis/source';
conditions = 'current_previous'; % 'previous' current_previous
if strcmp(conditions, 'current_previous')
    whichstat = 'statIC_ICvsC'; % statC_CvsN statIC_ICvsC statIC_ICvsN statC_CvsIC
elseif strcmp(conditions, 'previous')
    whichstat = 'statHemi';%'statHemi'; %statCvsN statICvsN
end

frequency = [10 20 40:4:100];
cnt = 0;
for k = frequency
  fprintf('computing T-statistic for frequency %d Hz\n', k);
  if strcmp(conditions, 'current_previous')
      d = dir(fullfile(datadir,sprintf('*3d4mm*pre_curprev_%d.mat',k)));
  elseif strcmp(conditions, 'previous')
      d = dir(fullfile(datadir,sprintf('*3d4mm*pre_%d.mat',k)));
  end
  for m = 1:numel(d)
    dum = load(fullfile(d(m).folder,d(m).name),whichstat);
    tmp(m) = dum.(whichstat);
  end
  clear dum
  dat = cat(2,tmp.stat);
  n   = size(dat,2);
  for m = 1:n
    dat(:,m+n) = nanmean(dat(:,m));
  end
  
  design   = [ones(1,n) ones(1,n)*2;1:n 1:n];
  cfg.ivar = 1;
  cfg.uvar = 2;
  tmp      = ft_statfun_depsamplesT(cfg, dat, design);
  cnt      = cnt+1;
  T(:,cnt) = tmp.stat;
  clear tmp;
  
end

load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
source = sourcemodel;
source.gamma1 = nanmean(T(:,3:8),2);
source.gamma2 = nanmean(T(:,10:16),2);
source.alpha  = T(:,1);
source.beta   = T(:,2);

cfgi = [];
cfgi.parameter = {'alpha' 'beta' 'gamma1' 'gamma2'};
source_int = ft_sourceinterpolate(cfgi, source, mri);

source_int.alpha(~isfinite(source_int.alpha))=0;
source_int.beta(~isfinite(source_int.beta))=0;
source_int.gamma1(~isfinite(source_int.gamma1))=0;
source_int.gamma2(~isfinite(source_int.gamma2))=0;

cmap = flipud(brewermap(64,'RdBu'));
cfgp=[];
cfgp.maskparameter = 'mask';
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'slice';
cfgp.nslices    =  30;
cfgp.slicerange =  [40 150];
cfgp.opacitylim = [2 4];

% cfgp.method = 'surface';

cfgp.funcolorlim = [-3 3];
cfgp.funparameter = 'alpha';
ft_sourceplot(cfgp, source_int);
title('alpha'); pause(0.01); 

cfgp.funcolorlim = [-3 3];
cfgp.funparameter = 'beta';
ft_sourceplot(cfgp, source_int);
title('beta'); pause(0.01); 

cfgp.funcolorlim = [-2 2];
cfgp.funparameter  = 'gamma1';
ft_sourceplot(cfgp, source_int);
title('gamma1'); pause(0.01); 
cfgp.funparameter  = 'gamma2';
ft_sourceplot(cfgp, source_int);
title('gamma2'); pause(0.01); 
