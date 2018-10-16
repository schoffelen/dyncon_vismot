function stat = vismot_spectral_prepost_stats(suffix, contrast, varargin)
% contrast: 10 item vector with (-)1's or zeros, depending on specific 
% contrast. items order as follows: [condition 1-5 pre, condition 1-5 post]
% condition 1: cue left, response left
% condition 2: cue left, response right
% condition 3: cue right, response left
% condition 4: cue right, response right
% condition 5: catch trial

dolog = ft_getopt(varargin , 'dolog', 0);

datadir = '/project/3011085.03/analysis/freq';
d   = dir(datadir);
sel = ~cellfun('isempty', strfind({d.name}', suffix));
d   = d(sel);

Fpre = cell(1,5);
Fpst = cell(1,5);
for k = 1:numel(d)
  fprintf('loading %s\n',fullfile(datadir,d(k).name));
  load(fullfile(datadir, d(k).name));
  for m = 1:5
    Fpre{k,m} = freqpre(m);
    Fpst{k,m} = freqpst(m);
  end
end

cfg           = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'powspctrm';
F = cell(1,10);
for k = 1:5
  F{k}   = ft_appendfreq(cfg, Fpre{:,k});
  F{k+5} = ft_appendfreq(cfg, Fpst{:,k});
end
clear Fpre Fpst;

if dolog
  cfg = [];
  cfg.operation = 'log10';
  cfg.parameter = 'powspctrm';
  for k = 1:numel(F)
    F{k} = ft_math(cfg, F{k});
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
      dat2.powspctrm = F{k}.powspctrm.*abs(contrast(k));
    else
      dat2.powspctrm = dat2.powspctrm + F{k}.powspctrm.*abs(contrast(k));
    end
    
    
  elseif sign(contrast(k))>0
    if ~exist('dat1', 'var')
      dat1 = F{k};
      dat1.powspctrm = F{k}.powspctrm.*abs(contrast(k));
    else
      dat1.powspctrm = dat1.powspctrm + F{k}.powspctrm.*abs(contrast(k));
    end
  end
end

n = size(dat1.powspctrm,1);
  
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 1000;
cfg.design = [ones(1,n) 2*ones(1,n);1:n 1:n];
stat = ft_freqstatistics(cfg, dat1, dat2);

