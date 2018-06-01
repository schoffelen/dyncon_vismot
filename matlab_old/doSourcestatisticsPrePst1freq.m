function doSourcestatisticsPrePst(band,doplot,nrand)

if nargin<3,
  nrand = 0;
end

if nargin<2,
  doplot = 0;
end

subjinfo;
%cd('/analyse/4/Project0030/source/4Dprepst');
cd('/analyse/1/Project0002/tmpProject0030/source/');

%exclude BKA01, GAR12 KBI24
names    = {SUBJ(:).name}';
selnames = find(~ismember(names,{'BKA01';'GAR12';'KBI24'}));
names    = names(selnames);

load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(names)
  %fname = [names{k}, 'stat4Dprepst',band];
  fname = [names{k}, 'stat',band,'prepst'];
  fprintf('loading %s\n', fname);
  load(fname);
  stat.pos = grid.pos;  %get positions consistent across subjects
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  tmp(stat.inside) = tmp(stat.inside) + 1;
  s{k,1}   = stat;
end
stat   = selectdata(s{:},  'param', {'stat2' 'stat2y' 'fwhm'});
nsubj  = size(stat.stat2,1);
inside = find(tmp==nsubj);
outside = setdiff(1:prod(dim), inside);
stat.inside  = inside(:);
stat.outside = outside(:);
clear tmp;

stat.stat2(:,outside,:)  = 0;
stat.stat2y(:,outside,:) = 0;

stat.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat.stat2y( nsubj+1:2*nsubj,:,:) = 0;
clear s;

stat.stat2   = double(stat.stat2);
stat.stat2y  = double(stat.stat2y);

%stat.stat2(nsubj+1:2*nsubj,inside,:)    = repmat(nanmean(stat.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat.stat2y(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.statistic = 'pooledT';
cfg.parameter = 'stat2y';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
%cfg.frequency = [4 36];
cfg.correctm  = 'max';
stat       = sourcestatistics(cfg, stat);
%cfg.frequency = [44 100];
%stathigh      = sourcestatistics(cfg, stat);

%cd('/analyse/4/Project0030/source/4Dprepst');
%save(['grandavg_',band], 'statlow', 'stathigh');
cd('/analyse/1/Project0002/tmpProject0030/source');
save(['grandavg_',band], 'stat');


if doplot,
  st13        = sourcestatistics(cfg, stat13);
  st42        = sourcestatistics(cfg, stat42);
  cfg.parameter = 'stat2y';
  st13y        = sourcestatistics(cfg, stat13);
  st42y        = sourcestatistics(cfg, stat42);
  cd /analyse/4/Project0030/figures/sourcedata/4Dpst
  for k = 1:length(names)
    fname = names{k};
    tmp1 = reshape(mean(stat42.stat2(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13.stat2(k,:,:),3),dim);
    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
    tmp1(outside) = nan;
    tmp1(1) = -max(abs(tmp1(:)));
    tmp1(2) =  max(abs(tmp1(:)));
    tmp2(outside) = nan;
    tmp2(1) = -max(abs(tmp2(:)));
    tmp2(2) =  max(abs(tmp2(:)));
    tmp(outside) = nan;
    tmp(1)  = -max(abs(tmp(:)));
    tmp(2)  =  max(abs(tmp(:)));
    figure;
    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
    cd Tstat
    print(gcf,'-dpng',[fname,band,'CongruencyStrat']);
    cd ..
    close
    tmp1 = reshape(mean(stat42.stat2y(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13.stat2y(k,:,:),3),dim);
    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
    tmp1(outside) = nan;
    tmp1(1) = -max(abs(tmp1(:)));
    tmp1(2) =  max(abs(tmp1(:)));
    tmp2(outside) = nan;
    tmp2(1) = -max(abs(tmp2(:)));
    tmp2(2) =  max(abs(tmp2(:)));
    tmp(outside) = nan;
    tmp(1)  = -max(abs(tmp(:)));
    tmp(2)  =  max(abs(tmp(:)));
    figure;
    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
    cd TstatYuen
    print(gcf,'-dpng',[fname,band,'CongruenctStrat']);
    cd ..
    close
  end
  
  %plot the averages
  tmp1 = reshape(st42.stat,dim);
  tmp2 = reshape(st13.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyTstat']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyRightTstat']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyLeftTstat']);
  close
  tmp1 = reshape(st42y.stat,dim);
  tmp2 = reshape(st13y.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyTstatYuen']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyRightTstatYuen']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyLeftTstatYuen']);
  close
end
%
%if nrand>0,
%  %addpath /home/jan/matlab/spm2;
%  s   = stat42;
%  s.stat2l = s.stat2;
%  s.stat2r = s.stat2;
%  s.stat2y = s.stat2;
%  s.stat2yl = s.stat2;
%  s.stat2yr = s.stat2;
%  dum = zeros(dim);
%  for k = 1:length(stat42.freq)
%    for m = 1:2*nsubj
%      dum(:) = reshape(stat42.stat2(m,:,k),dim);
%      spm_smooth(dum,dum,1.5);
%      s.stat2r(m,:,k) = dum(:);
%      dum(:) = reshape(stat13.stat2(m,:,k),dim);
%      spm_smooth(dum,dum,1.5);
%      s.stat2l(m,:,k) = dum(:);
%      dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
%      spm_smooth(dum,dum,1.5);
%      s.stat2(m,:,k) = dum(:);
%      dum(:) = (reshape(stat42.stat2y(m,:,k),dim)+flipdim(reshape(stat13.stat2y(m,:,k),dim),1))./sqrt(2);
%      spm_smooth(dum,dum,1.5);
%      s.stat2y(m,:,k) = dum(:);
%      dum(:) = reshape(stat42.stat2y(m,:,k),dim);
%      spm_smooth(dum,dum,1.5);
%      s.stat2yr(m,:,k) = dum(:);
%      dum(:) = reshape(stat13.stat2y(m,:,k),dim);
%      spm_smooth(dum,dum,1.5);
%      s.stat2yl(m,:,k) = dum(:);
%    end
%  end
%  sel = find(~ismember(names,{'GAR12' 'BKA01'}));
%  sel = [sel sel+nsubj];
%  cfg.design = cfg.design(:,sel);
%  cfg.design(2,:) = [1:numel(sel)/2 1:numel(sel)/2];
%  s   = selectdata(s,'rpt',sel);
%  cfg.numrandomization = nrand;
%  cfg.correctm         = 'max';
%  cfg.parameter        = 'stat2';
%  statC = sourcestatistics(cfg, s);
%  cfg.parameter = 'stat2l';
%  statCl = sourcestatistics(cfg, s);
%  cfg.parameter = 'stat2r';
%  statCr = sourcestatistics(cfg, s);
%  cd /analyse/4/Project0030/source/4Dpst
%  save(['grandavg_',band],'statC','statCl','statCr');
%  cfg.parameter        ='stat2y';
%  statC = sourcestatistics(cfg, s);
%  cfg.parameter = 'stat2yl';
%  statCl = sourcestatistics(cfg, s);
%  cfg.parameter = 'stat2yr';
%  statCr = sourcestatistics(cfg, s);
%  save(['grandavgYuen_',band],'statC','statCl','statCr');
%end
