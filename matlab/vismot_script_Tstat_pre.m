datadir = '/project/3011085.03/analysis/source';

frequency = [10 20 40:4:100];
cnt = 0;
for k = frequency
  fprintf('computing T-statistic for frequency %d Hz\n', k);
  d = dir(fullfile(datadir,sprintf('*3d4mm*pre*_%d.mat',k)));
  for m = 1:numel(d)
    load(fullfile(d(m).folder,d(m).name),'statResp');
    tmp(m) = statResp;
  end
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

cfgp.funparameter  = 'gamma1';
cfgp.maskparameter = 'mask';
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'slice';
cfgp.nslices    =  30;
cfgp.slicerange =  [40 150];
cfgp.opacitylim = [2 4];
cfgp.funcolorlim = [-4 4];

