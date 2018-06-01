function [stat,stat2,stat3] = collectVoxelstats(pathname,suffix,doplot,nrand,fnames)

if nargin<5,
  fnames = {'stat2' 'stat2y'};
end

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
dim = grid.dim;
for k = 1:length(names)
  fname = [names{k}, suffix];
  fprintf('loading %s\n', fname);
  load(fname);
  stat13.pos = grid.pos;  %get positions consistent across subjects
  stat42.pos = grid.pos;  %get positions consistent across subjects
  stat13.dim = grid.dim;
  stat42.dim = grid.dim;

  if k==1, dim = grid.dim; tmp = zeros(dim); end
  tmp(stat13.inside) = tmp(stat13.inside) + 1;
  tmp(stat42.inside) = tmp(stat42.inside) + 1;
  s{k,1}   = stat13;
  s2{k,1}  = stat42;
end
stat   = selectdata(s{:},   'param', fnames);
stat2  = selectdata(s2{:},  'param', fnames);
nsubj  = size(stat.(fnames{1}),1);
inside = find(tmp==nsubj*2);
outside = setdiff(1:prod(dim), inside);
stat.inside   = inside(:);
stat.outside  = outside(:);
stat2.inside  = inside(:);
stat2.outside = outside(:);
clear tmp;

for k = 1:numel(fnames)
  stat.(fnames{k})(:,outside,:) = 0;
  stat2.(fnames{k})(:,outside,:) = 0;
end
%stat.stat2(:,outside,:)   = 0;
%stat.stat2y(:,outside,:)  = 0;
%stat2.stat2(:,outside,:)  = 0;
%stat2.stat2y(:,outside,:) = 0;
stat3 = stat; 

fprintf('pooling the left/right data\n');
if ~isfield(stat, 'freq')
  stat.freq = 1:size(stat.(fnames{1}),3);
end
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    for j = 1:numel(fnames)
      tmp = (reshape(stat2.(fnames{j})(k,:,m), dim) + flipdim(reshape(stat.(fnames{j})(k,:,m), dim),1))./sqrt(2);
      stat3.(fnames{j})(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
    end
  end
end

fprintf('applying spatial smoothing\n');
smo = 1;
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    for j = 1:numel(fnames)
      tmp = reshape(stat.(fnames{j})(k,:,m), dim);
      spm_smooth(tmp, tmp, smo);
      stat.(fnames{j})(k,:,m) = tmp(:);
      tmp = reshape(stat2.(fnames{j})(k,:,m), dim);
      spm_smooth(tmp, tmp, smo);
      stat2.(fnames{j})(k,:,m) = tmp(:);
      tmp = reshape(stat3.(fnames{j})(k,:,m), dim);
      spm_smooth(tmp, tmp, smo);
      stat3.(fnames{j})(k,:,m) = tmp(:);
    end
  end
end  

for j = 1:numel(fnames)
  stat.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  stat.(fnames{j})   = double(stat.(fnames{j}));
  stat2.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  stat2.(fnames{j})   = double(stat2.(fnames{j}));
  stat3.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  stat3.(fnames{j})   = double(stat3.(fnames{j}));
end
clear s s2;

%stat.stat2(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean(stat.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat.stat2y(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat2.stat2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat2.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat2.stat2y(nsubj+1:2*nsubj,inside,:) = repmat(nanmean(stat2.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat3.stat2(nsubj+1:2*nsubj,inside,:)  = repmat(nanmean(stat3.stat2(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%stat3.stat2y(nsubj+1:2*nsubj,inside,:) = repmat(nanmean(stat3.stat2y(1:nsubj,inside,:),2),[1 numel(inside) 1]);

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = pos2dim3d(stat.pos);
cfg.inside    = stat.inside(:);
%cfg.statistic = 'pooledTtfce';
cfg.statistic = 'pooledT';
cfg.parameter = fnames{1};
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
stat3.dimord = 'rpt_pos_freq';
stat2.freq = stat.freq;
stat3.freq = stat.freq;
stat.dim = size(stat.(fnames{1}));;
stat2.dim = stat.dim;
stat3.dim = stat.dim;
statCl        = sourcestatistics(cfg, stat);
statCr        = sourcestatistics(cfg, stat2);
statC         = sourcestatistics(cfg, stat3);

cd(pathname);
save(['grandavg_',suffix,'_',fnames{1}],'statC','statCl','statCr');

%cfg.parameter        ='stat2y';
%statCl = sourcestatistics(cfg, stat);
%statCr = sourcestatistics(cfg, stat2);
%statC  = sourcestatistics(cfg, stat3);
%
%cd(pathname);
%save(['grandavgYuen_',suffix],'statC','statCl','statCr');
%
%if doplot,
%  st13        = sourcestatistics(cfg, stat13);
%  st42        = sourcestatistics(cfg, stat42);
%  cfg.parameter = 'stat2y';
%  st13y        = sourcestatistics(cfg, stat13);
%  st42y        = sourcestatistics(cfg, stat42);
%  cd /analyse/4/Project0030/figures/sourcedata/4Dpst
%  for k = 1:length(names)
%    fname = names{k};
%    tmp1 = reshape(mean(stat42.stat2(k,:,:),3),dim);
%    tmp2 = reshape(mean(stat13.stat2(k,:,:),3),dim);
%    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
%    tmp1(outside) = nan;
%    tmp1(1) = -max(abs(tmp1(:)));
%    tmp1(2) =  max(abs(tmp1(:)));
%    tmp2(outside) = nan;
%    tmp2(1) = -max(abs(tmp2(:)));
%    tmp2(2) =  max(abs(tmp2(:)));
%    tmp(outside) = nan;
%    tmp(1)  = -max(abs(tmp(:)));
%    tmp(2)  =  max(abs(tmp(:)));
%    figure;
%    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
%    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
%    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
%    cd Tstat
%    print(gcf,'-dpng',[fname,band,'CongruencyStrat']);
%    cd ..
%    close
%    tmp1 = reshape(mean(stat42.stat2y(k,:,:),3),dim);
%    tmp2 = reshape(mean(stat13.stat2y(k,:,:),3),dim);
%    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
%    tmp1(outside) = nan;
%    tmp1(1) = -max(abs(tmp1(:)));
%    tmp1(2) =  max(abs(tmp1(:)));
%    tmp2(outside) = nan;
%    tmp2(1) = -max(abs(tmp2(:)));
%    tmp2(2) =  max(abs(tmp2(:)));
%    tmp(outside) = nan;
%    tmp(1)  = -max(abs(tmp(:)));
%    tmp(2)  =  max(abs(tmp(:)));
%    figure;
%    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
%    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
%    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
%    cd TstatYuen
%    print(gcf,'-dpng',[fname,band,'CongruenctStrat']);
%    cd ..
%    close
%  end
%  
%  %plot the averages
%  tmp1 = reshape(st42.stat,dim);
%  tmp2 = reshape(st13.stat,dim);
%  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
%  tmp(1) = -max(abs(tmp(:)));
%  tmp(2) =  max(abs(tmp(:)));
%  tmp1(1) = -max(abs(tmp1(:)));
%  tmp1(2) =  max(abs(tmp1(:)));
%  tmp2(1) = -max(abs(tmp2(:)));
%  tmp2(2) =  max(abs(tmp2(:)));
%  figure;volplotJM(tmp,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyTstat']);
%  close
%  figure;volplotJM(tmp1,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyRightTstat']);
%  close
%  figure;volplotJM(tmp2,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyLeftTstat']);
%  close
%  tmp1 = reshape(st42y.stat,dim);
%  tmp2 = reshape(st13y.stat,dim);
%  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
%  tmp(1) = -max(abs(tmp(:)));
%  tmp(2) =  max(abs(tmp(:)));
%  tmp1(1) = -max(abs(tmp1(:)));
%  tmp1(2) =  max(abs(tmp1(:)));
%  tmp2(1) = -max(abs(tmp2(:)));
%  tmp2(2) =  max(abs(tmp2(:)));
%  figure;volplotJM(tmp,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyTstatYuen']);
%  close
%  figure;volplotJM(tmp1,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyRightTstatYuen']);
%  close
%  figure;volplotJM(tmp2,'montage');
%  print(gcf,'-dpng',[band,'AvgCongruencyLeftTstatYuen']);
%  close
%end
%
