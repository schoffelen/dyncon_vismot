function [allstat, stat] = collectVoxelstats2(pathname,suffix,doplot,nrand,fnames)

if nargin<5,
  fnames = {'stat' 'stat2'};
end

if ~iscell(fnames)
  fnames = {fnames};
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
%selnames = find(~ismember(names,{'PCL19';'BKA01';'GAR12';'KBI24'}));
names    = names(selnames);

load('/home/jan/matlab/mri/templategrid6mm.mat');
dim = grid.dim;
for k = 1:length(names)
  fname = [names{k}, suffix];
  fprintf('loading %s\n', fname);
  dat      = load(fname);
  varnames = fieldnames(dat);

  stat     = [];
  stat.pos = grid.pos;  %get positions consistent across subjects
  stat.dim = grid.dim;
  stat.inside  = double(dat.(varnames{1}).inside );
  stat.outside = double(dat.(varnames{1}).outside);

  newfnames = {};
  for m = 1:numel(varnames)
    for n = 1:numel(fnames)
      dum        = [varnames{m},fnames{n}];
      stat.(dum) = double(dat.(varnames{m}).(fnames{n}));
      stat.([dum,'dimord']) = 'pos';
      newfnames  = [newfnames {dum}];
    end
  end
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  tmp(stat.inside) = tmp(stat.inside) + 1;
  
  s{k,1} = stat;
  clear stat;
end

for k = 1:numel(newfnames)
  stat{k} = selectdata(s{:}, 'param', newfnames{k});
end
nsubj  = size(stat{1}.(newfnames{1}),1);
inside = find(tmp==nsubj);
outside = setdiff(1:prod(dim), inside);
for k = 1:numel(stat)
  stat{k}.inside   = inside(:);
  stat{k}.outside  = outside(:);
  stat{k}.dim = dim;
end
clear tmp;

for k = 1:numel(stat)
  stat{k}.(newfnames{k})(:,outside,:) = 0;
end

if ~isfield(stat{1}, 'freq')
  for k = 1:numel(stat) stat{k}.freq = 1:size(stat{k}.(newfnames{k}),3); end
end

fprintf('pooling the left/right data\n');
% try to establish which ones need to be pooled
for k = 1:numel(stat)
  for m = 1:numel(stat)
    tmp = (newfnames{k}~=newfnames{m});
    ndiff(k,m) = sum(tmp)==1 && (strcmp(newfnames{k}(tmp),'r') || strcmp(newfnames{k}(tmp),'l'));
  end
end
clear tmp

[ix,iy] = find(ndiff);
stat{5} = stat{1};
tmp     = fieldnames(stat{5});
for k = 1:numel(tmp)
  if size(stat{5}.(tmp{k}),1)==nsubj
    stat{5} = rmfield(stat{5},tmp{k});
  end
end
for k = 1:nsubj
  for m = 1:numel(stat{1}.freq)
    for j = 1:numel(newfnames)
      tmpname = [newfnames{ix(j)},'_',newfnames{iy(j)}];
      tmp = (reshape(stat{ix(j)}.(newfnames{ix(j)})(k,:,m), dim) + flipdim(reshape(stat{iy(j)}.(newfnames{iy(j)})(k,:,m), dim),1))./sqrt(2);
      stat{5}.(tmpname)(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
      stat{5}.([tmpname,'dimord']) = 'rpt_pos';
    end
  end
end

fprintf('applying spatial smoothing\n');
smo = 1.5;
for k = 1:numel(stat)
  tmp = fieldnames(stat{k});
  smoothfield = {};
  for m = 1:numel(tmp)
    if size(stat{k}.(tmp{m}),1)==nsubj
      smoothfield = [smoothfield tmp(m)];
    end
  end
  
  for i = 1:numel(smoothfield)
    tmpfield = stat{k}.(smoothfield{i});

    for m = 1:numel(stat{1}.freq)
      for j = 1:nsubj
        tmp = reshape(tmpfield(j,:,m), dim);
        spm_smooth(tmp, tmp, smo);
	tmpfield(j,:,m) = tmp(:);
      end
    end
    stat{k}.(smoothfield{i}) = tmpfield;
  end
end

fprintf('creating dummy volumes\n');
statfield = {};
statnum   = [];
for k = 1:numel(stat)
  tmp = fieldnames(stat{k});
  for m = 1:numel(tmp)
    if size(stat{k}.(tmp{m}),1)==nsubj
      stat{k}.(tmp{m})(nsubj+1:2*nsubj,:,:) = 0;
      %stat{k}.(tmp{m})(nsubj+1:2*nsubj,inside,:) = repmat(nanmedian(stat{k}.(tmp{m})(1:nsubj,inside,:), 2), [1 numel(inside) 1]);
      statfield = [statfield tmp(m)];
      statnum   = [statnum k];
    end
  end
end

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = dim; 
cfg.inside    = stat{1}.inside(:);
cfg.statistic = 'pooledT';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.correctm  = 'no';
cfg.implementation = 'new';
%cfg.correctm  = 'max';
if exist('frequency') && ~isempty(frequency)
  cfg.frequency = frequency;
else
  cfg.frequency = [stat{1}.freq(1) stat{1}.freq(end)];
end

for k = 1:numel(statnum)
  cfg.parameter = statfield{k};
  tmpstat = ft_sourcestatistics(cfg, stat{statnum(k)});
  allstat(k) = tmpstat;
end

cd(pathname);
save(['grandavg_',suffix],'allstat','stat');

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
