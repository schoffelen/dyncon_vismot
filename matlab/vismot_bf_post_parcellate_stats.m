function [stat, stat2] = vismot_bf_post_parcellate_stats(suffix, contrast, varargin)

dolog    = ft_getopt(varargin, 'dolog', 0);
fliphemi = ft_getopt(varargin, 'fliphemi', [0 0 0 0 0]);

datadir = '/home/language/jansch/projects/visuomotor/data/analyse/source';
d   = dir(datadir);
sel = ~cellfun('isempty', strfind({d.name}', suffix));
d   = d(sel);

load(fullfile('/home/language/jansch/projects/visuomotor/data/analyse/mri/Conte69_32k/atlas_subparc374_8k'));

cfg = [];
cfg.method = 'median';
cfg.parameter = 'powsmooth';

Fpst = cell(1,5);
for k = 1:numel(d)
  fprintf('loading %s\n',fullfile(datadir,d(k).name));
  load(fullfile(datadir, d(k).name));
  atlas.pos = source(1).pos;
	for m = 1:5
    Fpst{k,m} = ft_sourceparcellate(cfg, rmfield(source(m),{'trialinfo', 'cumtapcnt'}), atlas);
  end
end

for k = 1:numel(Fpst)
	Fpst{k}.dimord = 'chan_freq';
	Fpst{k}.pow = Fpst{k}.powsmooth;
end

cfg           = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'pow';
F = cell(1,5);
for k = 1:5
  F{k}   = ft_appendfreq(cfg, Fpst{:,k});
end
clear Fpst;

if dolog,
  cfg = [];
  cfg.operation = 'log10';
  cfg.parameter = 'pow';
  for k = 1:numel(F)
		if isfield(F{k}, 'tri'),
			tri  = F{k}.tri;
			F{k} = rmfield(F{k},'tri');
		end
    F{k} = ft_math(cfg, F{k});
		if exist('tri', 'var')
			F{k}.tri = tri;
		end
  end
end

if any(fliphemi),
  for k = 1:numel(fliphemi)
    if fliphemi(k)
      n = size(F{k}.label,1)./2;
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

n = size(dat1.pow,1);
  
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 2000;%5000;
cfg.correctm = 'no';
cfg.alpha = 0.025;
cfg.parameter = 'pow';
cfg.design = [ones(1,n) 2*ones(1,n);1:n 1:n];
stat = ft_freqstatistics(cfg, dat1, dat2);

stat.brainordinate = atlas;
stat2 = ft_checkdata(stat,'datatype','source');

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right