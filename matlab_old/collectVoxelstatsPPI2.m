function [stat,statV,statR] = collectVoxelstatsPPI2(pathname,suffix,doplot,nrand)

if nargin<4, nrand     = 0;            end
if nargin<3, doplot    = 0;            end
if nargin<2, error('need a suffix');   end
if nargin<1, error('need a pathname'); end

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
  stat.pos = grid.pos;  %get positions consistent across subjects
  stat.dim = [size(grid.pos,1) 1];
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  tmp(stat.inside) = tmp(stat.inside) + 1;
  s{k,1}   = stat;
end

fnames = fieldnames(s{1});
ok     = zeros(numel(fnames),1);
for k = 1:numel(fnames)
  tmp2 = size(s{1}.(fnames{k}));
  if isnumeric(s{1}.(fnames{k})) && numel(tmp2==2),
    ok(k,1) = all(tmp2==s{1}.dim);
  end
end
fnames = fnames(find(ok));

stat   = selectdata(s{:},   'param', fnames);
nsubj  = size(stat.stat2,1);
inside = find(tmp==nsubj);
outside = setdiff(1:prod(dim), inside);
stat.inside   = inside(:);
stat.outside  = outside(:);
clear tmp;

for k = 1:numel(fnames)
  stat.(fnames{k})(:,outside,:) = 0;
end

fieldpairs = {'roi'     'stat2'         'stat2b';
              'roiresp' 'stat2ppi1resp' 'stat2ppi2resp';
              'roivis'  'stat2ppi1vis'  'stat2ppi2vis'};
	       
fprintf('pooling the left/right data\n');
if ~isfield(stat, 'freq')
  stat.freq = 1:size(stat.stat2,3);
end
for j = 1:size(fieldpairs,1)
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    tmp = (reshape(stat.(fieldpairs{j,2})(k,:,m), dim) + flipdim(reshape(stat.(fieldpairs{j,3})(k,:,m), dim),1))./sqrt(2);
    stat.(fieldpairs{j,1})(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
  end
end
end

fnames = [fnames;{'roi';'roiresp';'roivis'}];

fprintf('applying spatial smoothing\n');
smo = 1.5;
for j = 1:numel(fnames)
for k = 1:nsubj
  for m = 1:numel(stat.freq)
    tmp = reshape(stat.(fnames{j})(k,:,m), dim);
    spm_smooth(tmp, tmp, smo);
    stat.(fnames{j})(k,:,m) = tmp(:);
  end
end  
end

for j = 1:numel(fnames)
  stat.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  stat.(fnames{j}) = double(stat.(fnames{j}));
end
clear s;

for j = 1:numel(fnames)
  stat.(fnames{j})(nsubj+1:2*nsubj,inside,:)   = repmat(nanmean(stat.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
end

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = pos2dim3d(stat.pos);
cfg.inside    = stat.inside(:);
%cfg.statistic = 'pooledT';
cfg.statistic = 'depsamplesT';
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
cfg.parameter = 'roivis';
stat.dimord = 'rpt_pos_freq';
stat.dim   = size(stat.(cfg.parameter));
statV      = sourcestatistics(cfg, stat);
cfg.parameter = 'roiresp';
statR      = sourcestatistics(cfg, stat);

cd(pathname);
save(['grandavg_',suffix],'stat','statV','statR');

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

