function [source, stat] = vismot_bf_pre(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency   = ft_getopt(varargin, 'frequency', 10);
smoothing   = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');
conditions  = ft_getopt(varargin, 'conditions', 'previous');
prewhiten   = ft_getopt(varargin, 'prewhiten', false);
lambda      = ft_getopt(varargin, 'lambda', '100%');

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
  end
end

freq =  vismot_spectral(subject,'output','fourier','toi', 'pre', 'conditions',conditions, 'foilim', [frequency frequency], 'smoothing', smoothing, 'prewhiten', prewhiten);
for k = 1:numel(freq)
  if ~isfield(freq(k),'trialinfo')
    freq(k).trialinfo = ones(numel(freq(k).cumtapcnt),1).*k;
  else
    if ~freq(k).trialinfo(:,end)==k
      freq(k).trialinfo(:,end+1) = k;
    end
  end
end

nconditions = numel(freq);
for k=1:nconditions
  tmpfreq{k} = freq(k);
end
cfg = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq = ft_appendfreq(cfg, tmpfreq{:});
clear tmpfreq
% freq = ft_appendfreq(cfg, freqpst(1), freqpst(2), freqpst(3), freqpst(4), freqpst(5));

% clear freqpst;

% load in the head model and the source model.
if isempty(sourcemodel)
  sourcemodel = vismot_anatomy_sourcemodel2d(subject);
end
load(fullfile(subject.pathname, 'headmodel', [subject.name, '_headmodel']));

% coregister the gradiometers if needed
if ~isempty(strfind(subject.datafile, 'h'))
  load(fullfile(subject.pathname, 'dewar2head_avg', [subject.name, 'dewar2head_avg']));
  
  % the transformation matrix M is in centimeters
  M(1:3,4)  = M(1:3,4)./100; % convert to meters
  freq.grad = ft_transform_geometry(M, ft_convert_units(freq.grad, 'm'));
end

headmodel   = ft_convert_units(headmodel,   'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');
if ~isfield(sourcemodel, 'inside')
  sourcemodel.inside = true(size(sourcemodel.pnt,1),1);
else
  try
    A = load('atlas_subparc374_8k.mat');
    idx = match_str(A.atlas.parcellationlabel, {'R_???_01', 'R_MEDIAL.WALL_01', 'L_???_01', 'L_MEDIAL.WALL_01'});
    sourcemodel.inside = sourcemodel.inside & ~ismember(A.atlas.parcellation, idx);
    %sourcemodel.inside(ismember(A.atlas.parcellation, idx))=0;
  end
end

% compute beamformer common spatial filters
cfg           = [];
cfg.grad      = freq.grad;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.channel   = freq.label;
cfg.singleshell.batchsize = 2000;
leadfield     = ft_prepare_leadfield(cfg);

[source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = sourcepow_pre(freq, headmodel, leadfield, lambda);

%%%%%% FIXME %%%%%%
%%%%%% then implement resampling %%%%%

function [source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = sourcepow_pre(freq, headmodel, leadfield, lambda)
if ~exist('onlycompute_stat13_42', 'var'); onlycompute_stat13_42=false; end

cfg                 = [];
cfg.grid            = leadfield;
cfg.headmodel       = headmodel;
cfg.method          = 'dics';
cfg.keeptrials      = 'yes';
cfg.dics.lambda     = lambda;
cfg.dics.fixedori   = 'yes';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfilter = 'yes';
tmpsource = ft_sourceanalysis(cfg, freq);
filter    = tmpsource.avg.filter;

cfg2             = [];
cfg2.fwhm        = 'yes';
if ~isfield(sourcemodel, 'dim')
  cfg2.fwhmmethod  = 'gaussfit';
  cfg2.fwhmmaxdist = 0.02;
end
fwhm             = ft_sourcedescriptives(cfg2, tmpsource);
fwhm             = fwhm.fwhm;

s     = keepfields(tmpsource,{'freq' 'tri' 'inside' 'pos' 'dim'});
if ~strcmp(conditions, 'current_previous')
  % same response hand contrast congruent minus incongruent
  stat13 = makesourcecontrast(freq, filter, s, [1 3], [], false, false, conditions);
  stat42 = makesourcecontrast(freq, filter, s, [4 2], [], false, false, conditions);
  
  if ~onlycompute_stat13_42
    statCvsIC  = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, true, conditions);
    statCvsIC2 = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, false, conditions);
    
    % same hemifield contrast congruent minus incongruent
    stat12 = makesourcecontrast(freq, filter, s, [1 2], [], false, false, conditions);
    stat43 = makesourcecontrast(freq, filter, s, [4 3], [], false, false, conditions);
    
    % vs neutral condition (don't stratify RT)
    stat15 = makesourcecontrast(freq, filter, s, [1 5], [], false, false, conditions);
    stat25 = makesourcecontrast(freq, filter, s, [2 5], [], false, false, conditions);
    stat35 = makesourcecontrast(freq, filter, s, [3 5], [], false, false, conditions);
    stat45 = makesourcecontrast(freq, filter, s, [4 5], [], false, false, conditions);
    
    % compute condition specific power
    cfg.grid.filter     = filter;
    cfg.dics.keepfilter = 'no';
    for k = 1:nconditions
      cfg2.trials = find(freq.trialinfo(:,end)==k);
      tmp         = ft_sourceanalysis(cfg, ft_selectdata(cfg2, freq));
      source(k)   = tmp;
    end
  end
else % not yet working FIXME
  % left response
  stat1_12 = makesourcecontrast(freq, filter, s, [1 1 2], [], false, false, conditions); % current C, previous C vs IC
  stat1_13 = makesourcecontrast(freq, filter, s, [1 1 3], [], false, false, conditions); % current C, previous C vs N
  stat3_45 = makesourcecontrast(freq, filter, s, [3 4 5], [], false, false, conditions); % current IC, previous IC vs C
  stat3_56 = makesourcecontrast(freq, filter, s, [3 5 6], [], false, false, conditions); % current IC, previous IC vs N
  % right response
  stat4_12 = makesourcecontrast(freq, filter, s, [4 1 2], [], false, false, conditions); % current C, previous C vs IC
  stat4_13 = makesourcecontrast(freq, filter, s, [4 1 3], [], false, false, conditions); % current C, previous C vs N
  stat2_45 = makesourcecontrast(freq, filter, s, [2 4 5], [], false, false, conditions); % current IC, previous IC vs C
  stat2_56 = makesourcecontrast(freq, filter, s, [2 5 6], [], false, false, conditions); % current IC, previous IC vs N
end


function stat = makesourcecontrast(freq, filter, s, contrast, whichflip, stratifyflag, poolhemi, conditions)

% compute contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest';

switch conditions
  case 'previous'
    % same response hand contrast
    % left hand response, C vs neutral
    tmpcfg              = [];
    tmpcfg.trials       = find(ismember(freq.trialinfo(:,end),[contrast whichflip]));
    tmpfreq           = ft_selectdata(tmpcfg, freq);
    s.pow = zeros(numel(s.inside),numel(tmpcfg.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    if ~isempty(whichflip)
      s = fliphemitrials(s, 'pow', 1, whichflip, contrast);
    end
    if poolhemi
      load standard_sourcemodel3d4mm; %FIXME this hardcoded assumes a 4mm sourcemodel!
      s = poolhemispheres(s, 'pow', 'left', [1 -1], sourcemodel);
    end
    
    if stratifyflag
      % stratify the trials based on the RTs
      RT = s.trialinfo(:,3);
      c  = s.trialinfo(:,1);
      
      tmpcfg        = [];
      tmpcfg.numbin = 5;
      out = ft_stratify(tmpcfg, RT(c==contrast(1))', RT(c==contrast(2))');
      s.trialinfo(c==contrast(1),3) = out{1}';
      s.trialinfo(c==contrast(2),3) = out{2}';
    end
    
    tmpcfg        = [];
    tmpcfg.trials = find(isfinite(s.trialinfo(:,3)));
    s             = ft_selectdata(tmpcfg, s);
    
    cfgs.design(s.trialinfo(:,4)==contrast(1)) = 1;
    cfgs.design(s.trialinfo(:,4)==contrast(2)) = 2;
    stat    = ft_sourcestatistics(cfgs, s);
    stat = rmfield(stat, {'prob', 'cirange', 'mask'});
    try, stat.tri = int16(stat.tri); end
    stat.pos = single(stat.pos);
    
  case 'current_previous'
    % left response
    % current C, previous C vs IC
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,1)==contrast(1) & freq.trialinfo(:,5),[contrast(2:3)]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    stat    = ft_sourcestatistics(cfgs, s);
    stat = rmfield(stat, {'prob', 'cirange', 'mask'});
    try, stat.tri = int16(stat.tri); end
    stat.pos = single(stat.pos);
end
