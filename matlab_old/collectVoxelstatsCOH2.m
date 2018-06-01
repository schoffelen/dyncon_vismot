function [statRi,statRc,statLi,statLc,roi1,roi2,statri,statrc,statli,statlc] = collectVoxelstatsCOH2(pathname,suffix,nrand,saveflag)

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
  
  statri.pos = grid.pos;  %get positions consistent across subjects
  statri.dim = [size(grid.pos,1) 1];
  statri.dcohdimord  = 'pos';  
  statri.dcohsdimord = 'pos';  
  statli.pos = grid.pos;  %get positions consistent across subjects
  statli.dim = [size(grid.pos,1) 1];
  statli.dcohdimord  = 'pos';  
  statli.dcohsdimord = 'pos';  
  statrc.pos = grid.pos;  %get positions consistent across subjects
  statrc.dim = [size(grid.pos,1) 1];
  statrc.dcohdimord  = 'pos';  
  statrc.dcohsdimord = 'pos';  
  statlc.pos = grid.pos;  %get positions consistent across subjects
  statlc.dim = [size(grid.pos,1) 1];
  statlc.dcohdimord  = 'pos';  
  statlc.dcohsdimord = 'pos';  
  
  roi1{k}  = statri.roi1;
  roi2{k}  = statri.roi2;
  
  if k==1, dim = grid.dim; tmp = zeros(dim); end
  if k==1, fwhm = zeros(dim); end 
  fwhm   = reshape(statri.fwhm, dim) + fwhm;
  statri = rmfield(statri, 'fwhm');
  tmp(statri.inside) = tmp(statri.inside) + 1;
  s{k,1}   = statri;
  s{k,2}   = statrc;
  s{k,3}   = statli;
  s{k,4}   = statlc;
  clear statri statrc statli statlc;
  rtr(k,:) = [s{k,1}.rt s{k,2}.rt];
  rtl(k,:) = [s{k,3}.rt s{k,4}.rt];
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
    statri = tmp; 
  else
    statri.(fnames{k}) = tmp.(fnames{k});
  end
  tmp = selectdata(s{:,2}, 'param', fnames{k});
  if k==1, 
    statrc = tmp;
  else
    statrc.(fnames{k}) = tmp.(fnames{k});
  end
  tmp = selectdata(s{:,3}, 'param', fnames{k});
  if k==1, 
    statli = tmp;
  else
    statli.(fnames{k}) = tmp.(fnames{k});
  end
  tmp = selectdata(s{:,4}, 'param', fnames{k});
  if k==1, 
    statlc = tmp;
  else
    statlc.(fnames{k}) = tmp.(fnames{k});
  end
end 

nsubj  = size(statri.dcohs,1);
inside = find(tmpinside==nsubj);
outside = setdiff(1:prod(dim), inside);
statri.inside   = inside(:);
statri.outside  = outside(:);
clear tmp;

for k = 1:numel(fnames)
  statri.(fnames{k})(:,outside,:) = 0;
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

%fprintf('applying spatial smoothing\n');
%smo = 2;
%for j = 1:numel(fnames)
%for k = 1:nsubj
%  for m = 1:numel(1) %FIXME works for 1 frequency only
%    tmp = reshape(statri.(fnames{j})(k,:,m), dim);
%    tmp(1) = tmp(2);
%    spm_smooth(tmp, tmp, smo);
%    statri.(fnames{j})(k,:,m) = tmp(:);
%    tmp = reshape(statrc.(fnames{j})(k,:,m), dim);
%    tmp(1) = tmp(2);
%    spm_smooth(tmp, tmp, smo);
%    statrc.(fnames{j})(k,:,m) = tmp(:);
%    tmp = reshape(statli.(fnames{j})(k,:,m), dim);
%    tmp(1) = tmp(2);
%    spm_smooth(tmp, tmp, smo);
%    statli.(fnames{j})(k,:,m) = tmp(:);
%    tmp = reshape(statlc.(fnames{j})(k,:,m), dim);
%    tmp(1) = tmp(2);
%    spm_smooth(tmp, tmp, smo);
%    statlc.(fnames{j})(k,:,m) = tmp(:);
%  end
%end  
%end

for j = 1:numel(fnames)
  statri.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statri.(fnames{j}) = double(statri.(fnames{j}));
  statrc.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statrc.(fnames{j}) = double(statrc.(fnames{j}));
  statli.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statli.(fnames{j}) = double(statli.(fnames{j}));
  statlc.(fnames{j})( nsubj+1:2*nsubj,:,:)   = 0;
  statlc.(fnames{j}) = double(statlc.(fnames{j}));
end
%clear s;

%for j = 1:numel(fnames)
%  statri.(fnames{j})(1:nsubj,inside,:)   = statri.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statri.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%  statrc.(fnames{j})(1:nsubj,inside,:)   = statrc.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statrc.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%  statli.(fnames{j})(1:nsubj,inside,:)   = statli.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statli.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%  statlc.(fnames{j})(1:nsubj,inside,:)   = statlc.(fnames{j})(1:nsubj,inside,:) - repmat(nanmean(statlc.(fnames{j})(1:nsubj,inside,:),2),[1 numel(inside) 1]);
%end

statri.dim = grid.dim;
statri.dcohsdimord  = statri.dcohdimord;
statrc.dim = grid.dim;
statrc.dcohsdimord  = statrc.dcohdimord;
statli.dim = grid.dim;
statli.dcohsdimord  = statli.dcohdimord;
statlc.dim = grid.dim;
statlc.dcohsdimord  = statlc.dcohdimord;

cfg           = [];
cfg.implementation = 'new';
cfg.method    = 'montecarlo';
cfg.numrandomization = nrand;
cfg.dim       = pos2dim3d(statri.pos);
cfg.inside    = statri.inside(:);
%cfg.statistic = 'pooledT';
cfg.statistic = 'yuenT';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.correctm  = 'max';

%cfg.parameter = 'dcohs';
cfg.parameter = 'dcoh';
statRi = ft_sourcestatistics(cfg, statri);
statRc = ft_sourcestatistics(cfg, statrc);
statLi = ft_sourcestatistics(cfg, statli);
statLc = ft_sourcestatistics(cfg, statlc);
%
%cfg.parameter = 'dcoh1';
%tmp = ft_sourcestatistics(cfg, stat);
%stat1.statu = tmp.stat;
%tmp = ft_sourcestatistics(cfg, statl);
%stat1l.statu = tmp.stat;
%tmp = ft_sourcestatistics(cfg, statr);
%stat1r.statu = tmp.stat;


statlc.rt = rtl;
statrc.rt = rtr;
statlc.fwhm = fwhm./17;
if saveflag,
  cd(pathname);
  save(['grandavg_',suffix],'statRi','statRc','statLc','statLi');
end

