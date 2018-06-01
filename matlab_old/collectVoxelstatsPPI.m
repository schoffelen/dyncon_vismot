function [stat,stat2,stat3] = collectVoxelstatsPPI(pathname,suffix,doplot,nrand)

if nargin<4,
  nrand     = 0;
end

if nargin<3,
  doplot = 0;
end

if nargin<2,
  error('need a suffix');
end

if nargin<1,
  error('need a pathname');
end

frequency = [];

subjinfo;
cd(pathname);

%exclude BKA01, GAR12 KBI24
names    = {SUBJ(:).name}';
selnames = find(~ismember(names,{'BKA01';'GAR12';'KBI24'}));
names    = names(selnames);

load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(names)
  fname = [names{k}, suffix];
  fprintf('loading %s\n', fname);
  load(fname);
  stat13.pos = grid.pos;  %get positions consistent across subjects
  stat12.pos = grid.pos;  %get positions consistent across subjects
  stat13.dim = [size(grid.pos,1) 1];
  stat12.dim = [size(grid.pos,1) 1];
  stat42.pos = grid.pos;  %get positions consistent across subjects
  stat43.pos = grid.pos;  %get positions consistent across subjects
  stat42.dim = [size(grid.pos,1) 1];
  stat43.dim = [size(grid.pos,1) 1];
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  tmp(stat13.inside) = tmp(stat13.inside) + 1;
  tmp(stat12.inside) = tmp(stat12.inside) + 1;
  tmp(stat42.inside) = tmp(stat42.inside) + 1;
  tmp(stat43.inside) = tmp(stat43.inside) + 1;
  s{k,1}   = stat13;
  s3{k,1}   = stat12;
  s2{k,1}  = stat42;
  s4{k,1}  = stat43;
end
stat   = selectdata(s{:},   'param', {'stat2' 'stat2b' 'stat2ppi1' 'stat2ppi2'});
stat2  = selectdata(s2{:},  'param', {'stat2' 'stat2b' 'stat2ppi1' 'stat2ppi2'});
stat3  = selectdata(s3{:},  'param', {'stat2' 'stat2b' 'stat2ppi1' 'stat2ppi2'});
stat4  = selectdata(s4{:},  'param', {'stat2' 'stat2b' 'stat2ppi1' 'stat2ppi2'});
nsubj  = size(stat.stat2,1);
inside = find(tmp==nsubj*4);
outside = setdiff(1:prod(dim), inside);
stat.inside   = inside(:);
stat.outside  = outside(:);
stat2.inside  = inside(:);
stat2.outside = outside(:);
stat3.inside  = inside(:);
stat3.outside = outside(:);
stat4.inside  = inside(:);
stat4.outside = outside(:);
clear tmp;

stat.stat2(:,outside,:)   = 0;
stat.stat2b(:,outside,:)  = 0;
stat2.stat2(:,outside,:)  = 0;
stat2.stat2b(:,outside,:) = 0;
stat3.stat2(:,outside,:)  = 0;
stat3.stat2b(:,outside,:) = 0;
stat4.stat2(:,outside,:)  = 0;
stat4.stat2b(:,outside,:) = 0;
stat.stat2ppi1(:,outside,:)  = 0;
stat.stat2ppi2(:,outside,:)  = 0;
stat2.stat2ppi1(:,outside,:) = 0;
stat2.stat2ppi2(:,outside,:) = 0;
stat3.stat2ppi1(:,outside,:) = 0;
stat3.stat2ppi2(:,outside,:) = 0;
stat4.stat2ppi1(:,outside,:) = 0;
stat4.stat2ppi2(:,outside,:) = 0;

statrsp = stat; 
statvis = stat3;

fprintf('pooling the left/right data\n');
if ~isfield(stat, 'freq')
  stat.freq = 1:size(stat.stat2,3);
end
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    tmp = (reshape(stat2.stat2(k,:,m), dim) + flipdim(reshape(stat.stat2b(k,:,m), dim),1))./sqrt(2);
    statrsp.stat2(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat2.stat2b(k,:,m), dim) + flipdim(reshape(stat.stat2(k,:,m), dim),1))./sqrt(2);
    statrsp.stat2b(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat2.stat2ppi1(k,:,m), dim) + flipdim(reshape(stat.stat2ppi2(k,:,m), dim),1))./sqrt(2);
    statrsp.stat2ppi1(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat2.stat2ppi2(k,:,m), dim) + flipdim(reshape(stat.stat2ppi1(k,:,m), dim),1))./sqrt(2);
    statrsp.stat2ppi2(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat4.stat2(k,:,m), dim) + flipdim(reshape(stat3.stat2b(k,:,m), dim),1))./sqrt(2);
    statvis.stat2(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat4.stat2b(k,:,m), dim) + flipdim(reshape(stat3.stat2(k,:,m), dim),1))./sqrt(2);
    statvis.stat2b(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat4.stat2ppi1(k,:,m), dim) + flipdim(reshape(stat3.stat2ppi2(k,:,m), dim),1))./sqrt(2);
    statvis.stat2ppi1(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    tmp = (reshape(stat4.stat2ppi2(k,:,m), dim) + flipdim(reshape(stat3.stat2ppi1(k,:,m), dim),1))./sqrt(2);
    statvis.stat2ppi2(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
  end
end

fprintf('applying spatial smoothing\n');
smo = 1;
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    tmp = reshape(stat.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat.stat2(k,:,m) = tmp(:);
    tmp = reshape(stat.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat.stat2b(k,:,m) = tmp(:);
    tmp = reshape(stat.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(stat.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat.stat2ppi2(k,:,m) = tmp(:);
    tmp = reshape(stat2.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat2.stat2(k,:,m) = tmp(:);
    tmp = reshape(stat2.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat2.stat2b(k,:,m) = tmp(:);
    tmp = reshape(stat2.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat2.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(stat2.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat2.stat2ppi2(k,:,m) = tmp(:);
    tmp = reshape(stat3.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat3.stat2(k,:,m) = tmp(:);
    tmp = reshape(stat3.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat3.stat2b(k,:,m) = tmp(:);
    tmp = reshape(stat3.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat3.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(stat4.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat4.stat2ppi2(k,:,m) = tmp(:);
    tmp = reshape(stat4.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat4.stat2(k,:,m) = tmp(:);
    tmp = reshape(stat4.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat4.stat2b(k,:,m) = tmp(:);
    tmp = reshape(stat4.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat4.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(stat4.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat4.stat2ppi2(k,:,m) = tmp(:);
    tmp = reshape(statrsp.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statrsp.stat2(k,:,m) = tmp(:);
    tmp = reshape(statrsp.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statrsp.stat2b(k,:,m) = tmp(:);
    tmp = reshape(statrsp.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statrsp.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(statrsp.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statrsp.stat2ppi2(k,:,m) = tmp(:);
    tmp = reshape(statvis.stat2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statvis.stat2(k,:,m) = tmp(:);
    tmp = reshape(statvis.stat2b(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statvis.stat2b(k,:,m) = tmp(:);
    tmp = reshape(statvis.stat2ppi1(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statvis.stat2ppi1(k,:,m) = tmp(:);
    tmp = reshape(statvis.stat2ppi2(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    statvis.stat2ppi2(k,:,m) = tmp(:);
  end
end  

stat.stat2( nsubj+1:2*nsubj,:,:)   = 0;
stat.stat2b( nsubj+1:2*nsubj,:,:)  = 0;
stat2.stat2( nsubj+1:2*nsubj,:,:)  = 0;
stat2.stat2b( nsubj+1:2*nsubj,:,:) = 0;
stat3.stat2( nsubj+1:2*nsubj,:,:)   = 0;
stat3.stat2b( nsubj+1:2*nsubj,:,:)  = 0;
stat4.stat2( nsubj+1:2*nsubj,:,:)  = 0;
stat4.stat2b( nsubj+1:2*nsubj,:,:) = 0;
statrsp.stat2( nsubj+1:2*nsubj,:,:)  = 0;
statrsp.stat2b( nsubj+1:2*nsubj,:,:) = 0;
statvis.stat2( nsubj+1:2*nsubj,:,:)  = 0;
statvis.stat2b( nsubj+1:2*nsubj,:,:) = 0;
stat.stat2ppi1( nsubj+1:2*nsubj,:,:)  = 0;
stat.stat2ppi2( nsubj+1:2*nsubj,:,:)  = 0;
stat2.stat2ppi1( nsubj+1:2*nsubj,:,:) = 0;
stat2.stat2ppi2( nsubj+1:2*nsubj,:,:) = 0;
stat3.stat2ppi1( nsubj+1:2*nsubj,:,:)  = 0;
stat3.stat2ppi2( nsubj+1:2*nsubj,:,:)  = 0;
stat4.stat2ppi1( nsubj+1:2*nsubj,:,:) = 0;
stat4.stat2ppi2( nsubj+1:2*nsubj,:,:) = 0;
statrsp.stat2ppi1( nsubj+1:2*nsubj,:,:) = 0;
statrsp.stat2ppi2( nsubj+1:2*nsubj,:,:) = 0;
statvis.stat2ppi1( nsubj+1:2*nsubj,:,:) = 0;
statvis.stat2ppi2( nsubj+1:2*nsubj,:,:) = 0;
clear s s2;

stat.stat2   = double(stat.stat2);
stat.stat2b  = double(stat.stat2b);
stat2.stat2  = double(stat2.stat2);
stat2.stat2b = double(stat2.stat2b);
stat3.stat2   = double(stat3.stat2);
stat3.stat2b  = double(stat3.stat2b);
stat4.stat2  = double(stat4.stat2);
stat4.stat2b = double(stat4.stat2b);
statrsp.stat2  = double(statrsp.stat2);
statrsp.stat2b = double(statrsp.stat2b);
statvis.stat2  = double(statvis.stat2);
statvis.stat2b = double(statvis.stat2b);
stat.stat2ppi1   = double(stat.stat2ppi1);
stat.stat2ppi2   = double(stat.stat2ppi2);
stat2.stat2ppi1  = double(stat2.stat2ppi1);
stat2.stat2ppi2  = double(stat2.stat2ppi2);
stat3.stat2ppi1   = double(stat3.stat2ppi1);
stat3.stat2ppi2   = double(stat3.stat2ppi2);
stat4.stat2ppi1  = double(stat4.stat2ppi1);
stat4.stat2ppi2  = double(stat4.stat2ppi2);
statrsp.stat2ppi1  = double(statrsp.stat2ppi1);
statrsp.stat2ppi2  = double(statrsp.stat2ppi2);
statvis.stat2ppi1  = double(statvis.stat2ppi1);
statvis.stat2ppi2  = double(statvis.stat2ppi2);

stat.stat2ppi1(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean( stat.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat.stat2ppi2(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean( stat.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat2.stat2ppi1(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat2.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat2.stat2ppi2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat2.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat3.stat2ppi1(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean( stat3.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat3.stat2ppi2(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean( stat3.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat4.stat2ppi1(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat4.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
stat4.stat2ppi2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat4.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
statrsp.stat2ppi1(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(statrsp.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
statrsp.stat2ppi2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(statrsp.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
statvis.stat2ppi1(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(statvis.stat2ppi1(1:nsubj,inside,:),2),[1 numel(inside) 1]);
statvis.stat2ppi2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(statvis.stat2ppi2(1:nsubj,inside,:),2),[1 numel(inside) 1]);

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = pos2dim3d(stat.pos);
cfg.inside    = stat.inside(:);
%cfg.statistic = 'pooledTtfce';
cfg.statistic = 'pooledT';
cfg.parameter = 'stat2ppi1';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.correctm  = 'no';
%cfg.correctm  = 'max';
%cfg.correctm  = 'cluster';
%cfg.clusterthreshold = 'nonparametric_individual';
%cfg.clusteralpha = 0.025;
%cfg.dim = pos2dim3d(stat.pos);
%cfg.clustercritval = [-2 2];
%cfg.clusterthreshold = 'parametric';
if ~isempty(frequency)
  cfg.frequency = frequency;
else
  cfg.frequency = [stat.freq(1) stat.freq(end)];
end
stat.dimord = 'rpt_pos_freq';
stat2.dimord = 'rpt_pos_freq';
statrsp.dimord = 'rpt_pos_freq';
statvis.dimord = 'rpt_pos_freq';
stat2.freq = stat.freq;
statrsp.freq = stat.freq;
statvis.freq = stat.freq;
stat.dim   = size(stat.(cfg.parameter));
stat2.dim   = size(stat.(cfg.parameter));
stat3.dim   = size(stat.(cfg.parameter));
stat4.dim   = size(stat.(cfg.parameter));
statrsp.dim   = size(stat.(cfg.parameter));
statvis.dim   = size(stat.(cfg.parameter));
statCl        = sourcestatistics(cfg, stat);
statCr        = sourcestatistics(cfg, stat2);
statC         = sourcestatistics(cfg, statrsp);
statVl        = sourcestatistics(cfg, stat3);
statVr        = sourcestatistics(cfg, stat4);
statV         = sourcestatistics(cfg, statvis); %FIXME also do stat3 stat4

cfg.parameter = 'stat2ppi2';
statCl2        = sourcestatistics(cfg, stat);
statCr2        = sourcestatistics(cfg, stat2);
statC2         = sourcestatistics(cfg, statrsp);
statVl2        = sourcestatistics(cfg, stat3);
statVr2        = sourcestatistics(cfg, stat4);
statV2         = sourcestatistics(cfg, statvis); %FIXME also do stat3 stat4

cd(pathname);
save(['grandavg_',suffix],'statC','statCl','statCr','statC2','statCl2','statCr2','statV','statVl','statVr','statV2','statVl2','statVr2','stat','stat2','stat3','stat4','statrsp','statvis');

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

