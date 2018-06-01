function [stat,statl,statr,roi1,roi2,stat1l,stat2l,stat1r,stat2r,stat1,stat2] = collectVoxelstatsCOH(pathname,suffix,nrand,saveflag)

if nargin<4, saveflag  = 0;            end
if nargin<3, nrand     = 0;            end
if nargin<2, error('need a suffix');   end
if nargin<1, error('need a pathname'); end

frequency = [];

subjinfo;
cd(pathname);

%exclude BKA01, GAR12 KBI24
names    = {SUBJ(:).name}';
%selnames = find(~ismember(names,{'BKA01';'GAR12';'KBI24'}));
selnames = find(~ismember(names,{'PCL19';'VIA12';'BKA01';'GAR12';'KBI24'}));
names    = names(selnames);

load('/home/jan/matlab/mri/templategrid6mm.mat');
%names = names([1 2 4 5 7 8 10:numel(names)]);
%names = names([
for k = 1:length(names)
  fname = [names{k}, suffix];
  fprintf('loading %s\n', fname);
  load(fname);
  
  stat.pos = grid.pos;  %get positions consistent across subjects
  stat.dim = [size(grid.pos,1) 1];
  stat.dcoh1dimord = 'pos';  
  stat.dcoh2dimord = 'pos';  
  stat.dcoh1sdimord = 'pos';  
  stat.dcoh2sdimord = 'pos';  
  
  roi1{k}  = stat.roi1;
  roi2{k}  = stat.roi2;
  
  statl.pos = grid.pos;  %get positions consistent across subjects
  statl.dim = [size(grid.pos,1) 1];
  statl.dcoh1dimord = 'pos';  
  statl.dcoh2dimord = 'pos';  
  statl.dcoh1sdimord = 'pos';  
  statl.dcoh2sdimord = 'pos';  

  statr.pos = grid.pos;  %get positions consistent across subjects
  statr.dim = [size(grid.pos,1) 1];
  statr.dcoh1dimord = 'pos';  
  statr.dcoh2dimord = 'pos';  
  statr.dcoh1sdimord = 'pos';  
  statr.dcoh2sdimord = 'pos';  
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  if k==1, fwhm = zeros(dim); end 
  fwhm = reshape(stat.fwhm, dim) + fwhm;
  stat = rmfield(stat, 'fwhm');
  tmp(stat.inside) = tmp(stat.inside) + 1;
  s{k,1}   = stat;
  s{k,2}   = statl;
  s{k,3}   = statr;
  clear stat statl statr;
  rtl(k,:) = s{k,2}.rt;
  rtr(k,:) = s{k,3}.rt;
end
tmpinside = tmp;

fnames = fieldnames(s{1});
ok     = zeros(numel(fnames),1);
for k = 1:numel(fnames)
  tmp2 = size(s{1}.(fnames{k}));
  if isnumeric(s{1}.(fnames{k})) && numel(tmp2==2),
    ok(k,1) = all(tmp2==s{1}.dim);
  end
end
fnames = fnames(find(ok));

for k = 1:numel(fnames)
  tmp = selectdata(s{:,1}, 'param', fnames{k});
  if k==1, 
    stat = tmp; 
  else
    stat.(fnames{k}) = tmp.(fnames{k});
  end
  tmp = selectdata(s{:,2}, 'param', fnames{k});
  if k==1, 
    statl = tmp;
  else
    statl.(fnames{k}) = tmp.(fnames{k});
  end
  tmp = selectdata(s{:,3}, 'param', fnames{k});
  if k==1, 
    statr = tmp;
  else
    statr.(fnames{k}) = tmp.(fnames{k});
  end
end 

nsubj  = size(stat.dcoh1s,1);
inside = find(tmpinside==nsubj);
outside = setdiff(1:prod(dim), inside);
stat.inside   = inside(:);
stat.outside  = outside(:);
clear tmp;

for k = 1:numel(fnames)
  stat.(fnames{k})(:,outside,:) = 0;
end

%fieldpairs = {'roi'     'stat2'         'stat2b';
%              'roiresp' 'stat2ppi1resp' 'stat2ppi2resp';
%              'roivis'  'stat2ppi1vis'  'stat2ppi2vis'};
%	       
%fprintf('pooling the left/right data\n');
%if ~isfield(stat, 'freq')
%  stat.freq = 1:size(stat.stat2,3);
%end
%for j = 1:size(fieldpairs,1)
%for k = 1:nsubj
%  for m = 1:numel(stat.freq)
%    tmp = (reshape(stat.(fieldpairs{j,2})(k,:,m), dim) + flipdim(reshape(stat.(fieldpairs{j,3})(k,:,m), dim),1))./sqrt(2);
%    stat.(fieldpairs{j,1})(k,:,m) = reshape(tmp, [1 numel(tmp) 1]);
%  end
%end
%end

fprintf('applying spatial smoothing\n');
smo = 2;
for j = 1:numel(fnames)
for k = 1:nsubj
  for m = 1:numel(1) %FIXME works for 1 frequency only
    tmp = reshape(stat.(fnames{j})(k,:,m), dim);
    tmp(1) = tmp(2);
    spm_smooth(tmp, tmp, smo);
    stat.(fnames{j})(k,:,m) = tmp(:);
    tmp = reshape(statl.(fnames{j})(k,:,m), dim);
    tmp(1) = tmp(2);
    spm_smooth(tmp, tmp, smo);
    statl.(fnames{j})(k,:,m) = tmp(:);
    tmp = reshape(statr.(fnames{j})(k,:,m), dim);
    tmp(1) = tmp(2);
    spm_smooth(tmp, tmp, smo);
    statr.(fnames{j})(k,:,m) = tmp(:);
  end
end  
end

for j = 1:numel(fnames)
  stat.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  stat.(fnames{j}) = double(stat.(fnames{j}));
  statl.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statl.(fnames{j}) = double(statl.(fnames{j}));
  statr.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statr.(fnames{j}) = double(statr.(fnames{j}));
end
%clear s;

for j = 1:numel(fnames)
  stat.(fnames{j})(1:nsubj,inside,:)   = stat.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(stat.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
  statl.(fnames{j})(1:nsubj,inside,:)   = statl.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statl.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
  statr.(fnames{j})(1:nsubj,inside,:)   = statr.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statr.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
end

stat.dim = grid.dim;
statl.dim = grid.dim;
statr.dim = grid.dim;
stat.dcoh1dimord  = stat.dcoh1dimord;
stat.dcoh2dimord  = stat.dcoh1dimord;
statl.dcoh1dimord = stat.dcoh1dimord;
statl.dcoh2dimord = stat.dcoh1dimord;
statr.dcoh1dimord = stat.dcoh1dimord;
statr.dcoh2dimord = stat.dcoh1dimord;
stat.dcoh1sdimord  = stat.dcoh1dimord;
stat.dcoh2sdimord  = stat.dcoh1dimord;
statl.dcoh1sdimord = stat.dcoh1dimord;
statl.dcoh2sdimord = stat.dcoh1dimord;
statr.dcoh1sdimord = stat.dcoh1dimord;
statr.dcoh2sdimord = stat.dcoh1dimord;

cfg           = [];
cfg.implementation = 'new';
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = pos2dim3d(stat.pos);
cfg.inside    = stat.inside(:);
%cfg.statistic = 'pooledT';
cfg.statistic = 'yuenT';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.correctm  = 'max';

cfg.parameter = 'dcoh1s';
stat1 = ft_sourcestatistics(cfg, stat);
%stat1 = [];
stat1l = ft_sourcestatistics(cfg, statl);
stat1r = ft_sourcestatistics(cfg, statr);

cfg.parameter = 'dcoh1';
tmp = ft_sourcestatistics(cfg, stat);
stat1.statu = tmp.stat;
tmp = ft_sourcestatistics(cfg, statl);
stat1l.statu = tmp.stat;
tmp = ft_sourcestatistics(cfg, statr);
stat1r.statu = tmp.stat;

cfg.parameter = 'dcoh2s';
stat2 = ft_sourcestatistics(cfg, stat);
%stat2 = [];
stat2l = ft_sourcestatistics(cfg, statl);
stat2r = ft_sourcestatistics(cfg, statr);

cfg.parameter = 'dcoh2';
tmp = ft_sourcestatistics(cfg, stat);
stat2.statu = tmp.stat;
tmp = ft_sourcestatistics(cfg, statl);
stat2l.statu = tmp.stat;
tmp = ft_sourcestatistics(cfg, statr);
stat2r.statu = tmp.stat;

statl.rt = rtl;
statr.rt = rtr;
stat.fwhm = fwhm./17;
if saveflag,
  cd(pathname);
  save(['grandavg_',suffix],'stat1','stat1l','stat1r','stat2','stat2l','stat2r');
end

