function doSourcestatisticsPst(band,doplot,nrand)

if nargin<3,
  nrand = 0;
end

if nargin<2,
  doplot = 0;
end

cd('/analyse/4/Project0030/source/4Dpst');

d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,band)) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
    fprintf('loading %s\n', fname);
    load(fname);
    stat13.pos  = grid.pos;
    stat42.pos  = grid.pos;
    cnt       = cnt+1;
    s13{cnt}  = stat13;
    s42{cnt}  = stat42;
    names{end+1} = fname(1:5);
  end
end
stat13  = selectdata(s13{:},  'param', {'stat2' 'stat2y'});
stat42  = selectdata(s42{:},  'param', {'stat2' 'stat2y'});
nsubj   = size(stat13.stat2,1);
stat13.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat42.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat13.stat2y( nsubj+1:2*nsubj,:,:) = 0;
stat42.stat2y( nsubj+1:2*nsubj,:,:) = 0;
clear s13 s42;

tmp = stat13;
tmp = fixinside(tmp, 'logical');
dim = pos2dim3d(tmp.pos);
inside = reshape(tmp.inside, dim) + flipdim(reshape(tmp.inside, dim),1);
tmp.inside = inside==2;
tmp = fixinside(tmp, 'index');
inside  = tmp.inside;
outside = tmp.outside;
clear tmp;

stat13.stat2(:,outside,:) = 0;
stat13.stat2y(:,outside,:) = 0;
stat42.stat2(:,outside,:) = 0;
stat42.stat2y(:,outside,:) = 0;
stat13.inside = inside;
stat13.outside = outside;
stat42.inside = inside;
stat42.outside = outside;

stat13.stat2   = double(stat13.stat2);
stat13.stat2y  = double(stat13.stat2y);
stat42.stat2   = double(stat42.stat2);
stat42.stat2y  = double(stat42.stat2y);

stat13.stat2(nsubj+1:2*nsubj,inside,:)    = repmat(nanmean(stat13.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat13.stat2y(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat13.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat42.stat2(nsubj+1:2*nsubj,inside,:)    = repmat(nanmean(stat42.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat42.stat2y(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat42.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'pooledT';
%cfg.statistic = 'depsamplesT';
cfg.parameter = 'stat2';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.frequency = [stat13.freq(1) stat13.freq(end)];
cfg.avgoverfreq = 'yes';

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

if nrand>0,
  %addpath /home/jan/matlab/spm2;
  s   = stat42;
  s.stat2l = s.stat2;
  s.stat2r = s.stat2;
  s.stat2y = s.stat2;
  s.stat2yl = s.stat2;
  s.stat2yr = s.stat2;
  dum = zeros(dim);
  for k = 1:length(stat42.freq)
    for m = 1:2*nsubj
      dum(:) = reshape(stat42.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2r(m,:,k) = dum(:);
      dum(:) = reshape(stat13.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2l(m,:,k) = dum(:);
      dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2(m,:,k) = dum(:);
      dum(:) = (reshape(stat42.stat2y(m,:,k),dim)+flipdim(reshape(stat13.stat2y(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2y(m,:,k) = dum(:);
      dum(:) = reshape(stat42.stat2y(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2yr(m,:,k) = dum(:);
      dum(:) = reshape(stat13.stat2y(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2yl(m,:,k) = dum(:);
    end
  end
  sel = find(~ismember(names,{'GAR12' 'BKA01'}));
  sel = [sel sel+nsubj];
  cfg.design = cfg.design(:,sel);
  cfg.design(2,:) = [1:numel(sel)/2 1:numel(sel)/2];
  s   = selectdata(s,'rpt',sel);
  cfg.numrandomization = nrand;
  cfg.correctm         = 'max';
  cfg.parameter        = 'stat2';
  statC = sourcestatistics(cfg, s);
  cfg.parameter = 'stat2l';
  statCl = sourcestatistics(cfg, s);
  cfg.parameter = 'stat2r';
  statCr = sourcestatistics(cfg, s);
  cd /analyse/4/Project0030/source/4Dpst
  save(['grandavg_',band],'statC','statCl','statCr');
  cfg.parameter        ='stat2y';
  statC = sourcestatistics(cfg, s);
  cfg.parameter = 'stat2yl';
  statCl = sourcestatistics(cfg, s);
  cfg.parameter = 'stat2yr';
  statCr = sourcestatistics(cfg, s);
  save(['grandavgYuen_',band],'statC','statCl','statCr');
end
