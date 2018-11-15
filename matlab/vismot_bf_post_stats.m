function stat = vismot_bf_post_stats(suffix, contrast, varargin)

dolog    = ft_getopt(varargin, 'dolog', 0);
fliphemi = ft_getopt(varargin, 'fliphemi', [0 0 0 0 0]);

datadir = '/project/3011085.03/analysis/source/mve';
d   = dir(datadir);
sel = ~cellfun('isempty', strfind({d.name}', suffix));
d   = d(sel);

Fpst = cell(1,5);
for k = 1:numel(d)
  fprintf('loading %s\n',fullfile(datadir,d(k).name));
  load(fullfile(datadir, d(k).name));
  for m = 1:5
    Fpst{k,m} = source(m);
  end
end

% sweep across all data and explicitly impose the same pos-field everywhere
% (this is allowed on the assumption that the data is coregistered in some
% normalized space
for k = 1:numel(Fpst)
  Fpst{k}.pos = Fpst{1}.pos;
end

cfg           = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'powsmooth';
F = cell(1,5);
for k = 1:5
  F{k}   = ft_appendsource(cfg, Fpst{:,k});
	F{k}.pow = F{k}.powsmooth;
end
clear Fpst;

if dolog
  cfg = [];
  cfg.operation = 'log10';
  cfg.parameter = 'pow';
  for k = 1:numel(F)
		if isfield(F{k}, 'tri')
			tri  = F{k}.tri;
			F{k} = rmfield(F{k},'tri');
		end
    F{k} = ft_math(cfg, F{k});
		if exist('tri', 'var')
			F{k}.tri = tri;
		end
  end
end

if any(fliphemi)
  for k = 1:numel(fliphemi)
    if fliphemi(k)
      n = size(F{k}.pos,1)./2;
      F{k}.pow = F{k}.pow(:,[n+(1:n) 1:n],:);
    end
  end
end

contrast(contrast<0) = contrast(contrast<0)./sum(abs(contrast(contrast<0)));
contrast(contrast>0) = contrast(contrast>0)./sum(contrast(contrast>0));

for k = 1:numel(contrast)
  if contrast(k)==0,
    % don't use this data
  elseif sign(contrast(k))<0
    if ~exist('dat2', 'var')
      dat2 = F{k};
      dat2.pow = F{k}.pow.*abs(contrast(k));
    else
      dat2.pow = dat2.pow + F{k}.pow.*abs(contrast(k));
    end
    
    
  elseif sign(contrast(k))>0
    if ~exist('dat1', 'var')
      dat1 = F{k};
      dat1.pow = F{k}.pow.*abs(contrast(k));
    else
      dat1.pow = dat1.pow + F{k}.pow.*abs(contrast(k));
    end
  end
end
dat1.inside = isfinite(dat1.pow(1,:))';
dat2.inside = dat1.inside;

n = size(dat1.pow,1);
  
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 2000;%5000;
cfg.correctm = 'no';
cfg.alpha = 0.025;
cfg.design = [ones(1,n) 2*ones(1,n);1:n 1:n];
stat = ft_sourcestatistics(cfg, dat1, dat2);

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right