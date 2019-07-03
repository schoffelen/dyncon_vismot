
function [source, stat, filter] = vismot_bf(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency    = ft_getopt(varargin, 'frequency', 10);
smoothing    = ft_getopt(varargin, 'smoothing', []);
sourcemodel  = ft_getopt(varargin, 'sourcemodel');
prewhiten    = istrue(ft_getopt(varargin, 'prewhiten', false));
lambda       = ft_getopt(varargin, 'lambda', ' 100%');
nrand        = ft_getopt(varargin, 'nrand', 0); % number of randomization for sensor subsampling
N            = ft_getopt(varargin, 'N', 90);
latoi        = ft_getopt(varargin, 'latoi', []);
conditions   = ft_getopt(varargin, 'conditions', []);
stratifyflag = ft_getopt(varargin, 'stratifyflag', false);

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
  end
end
if isempty(latoi)
  error('please specify latoi (pre or post)')
end

% post/pre specific settings
if strcmp(latoi, 'post')
  toi = 'prepost';
  if isempty(conditions); conditions = 'current'; end
elseif strcmp(latoi, 'pre')
  toi = 'pre';
  if isempty(conditions); conditions = 'previous'; end
end

freq    = vismot_spectral(subject, 'foilim', [frequency frequency], 'toi', toi, 'balance', true, 'smoothing', smoothing, 'output', 'fourier', 'prewhiten', prewhiten);

if strcmp(latoi, 'post')
  for k = 1:5
    % add marker for pre/pst
    freq(1,k).trialinfo(:,end+1) = 1;
    freq(2,k).trialinfo(:,end+1) = 2;
  end
elseif strcmp(latoi, 'pre')
  for k = 1:numel(freq)
    if ~isfield(freq(k),'trialinfo')
      freq(k).trialinfo = ones(numel(freq(k).cumtapcnt),1).*k;
    else
      if ~freq(k).trialinfo(:,end)==k
        freq(k).trialinfo(:,end+1) = k;
      end
    end
  end
end

for k=1:size(freq,2)
  tmpfreq{k} = freq(1,k);
  try
    tmpfreq{k+size(freq,2)} = freq(2,k);
  catch
  end
end
cfg = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq = ft_appendfreq(cfg, tmpfreq{:});
clear tmpfreq

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

% compute leadfields and beamformer common spatial filters
cfg           = [];
cfg.grad      = freq.grad;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.channel   = freq.label;
cfg.singleshell.batchsize = 2000;
leadfield     = ft_prepare_leadfield(cfg);

[source, stat, filter] = sourcepow(freq, headmodel, sourcemodel, leadfield, lambda, latoi, conditions, stratifyflag);

if nrand>0
  
  dpow13 = single(zeros(size(stat.stat13.stat)));
  dpow42 = dpow13;
  dpow13sq = dpow13;
  dpow42sq = dpow42;
  
  for k=1:nrand
    tic;
    indx = sort(randperm(numel(freq(1).label),N)); % keep it sorted!! -> subsampling of sensors
    tmpleadfield = leadfield;
    tmpleadfield.leadfield(tmpleadfield.inside) = cellrowselect(tmpleadfield.leadfield(tmpleadfield.inside),indx);
    
    tmpcfg = [];
    tmpcfg.channel = freq(1).label(indx);
    for m = 1:numel(freq)
      tmpfreq(m) = ft_selectdata(tmpcfg, freq(m));
    end
    tmpleadfield.label = tmpfreq(1).label;
    
    tmpsource = sourcepow(tmpfreq, headmodel, sourcemodel, tmpleadfield, lambda, latoi, conditions, stratifyflag, true);
    
    dpow13   = dpow13   + double(tmpsource(1).avg.pow - tmpsource(3).avg.pow);
    dpow13sq = dpow13sq + double(tmpsource(1).avg.pow - tmpsource(3).avg.pow).^2;
    dpow42   = dpow42   + double(tmpsource(4).avg.pow - tmpsource(2).avg.pow);
    dpow42sq = dpow42sq + double(tmpsource(4).avg.pow - tmpsource(2).avg.pow).^2;
    looptime(k) = toc;
  end
  
  
  mx13 = dpow13./nrand;
  sx13 = sqrt((dpow13sq-(dpow13.^2)./nrand)./nrand);
  zx13 = mx13./sx13;
  
  mx42 = dpow42./nrand;
  sx42 = sqrt((dpow42sq-(dpow42.^2)./nrand)./nrand);
  zx42 = mx42./sx42;
  
  stat.zx13 = zx13;
  stat.zx42 = zx42;
end


function [source, stat, filter] = sourcepow(freq, headmodel, sourcemodel, leadfield, lambda, latoi, conditions, stratifyflag,  onlycompute_stat13_42)
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

s     = keepfields(tmpsource, {'freq' 'tri' 'inside' 'pos' 'dim'});

% compute condition specific power, this is without stratification for RT
cfg.grid.filter     = filter;
cfg.dics.keepfilter = 'no';
for k = 1:5
  if strcmp(latoi, 'post')
    tmpcfg.trials = find(freq.trialinfo(:,1)==k & freq.trialinfo(:,end)==2); % for the pst trials only
  elseif strcmp(latoi, 'pre')
    tmpcfg.trials = find(freq.trialinfo(:,end)==k);
  end
  source(k)     = ft_sourceanalysis(cfg, ft_selectdata(tmpcfg, freq));
end

if ~onlycompute_stat13_42
  if ~strcmp(conditions, 'current_previous')
    % same response hand contrast congruent minus incongruent
    stat13 = makesourcecontrast(freq, filter, s, [1 3], [], stratifyflag, false, latoi, conditions);
    stat42 = makesourcecontrast(freq, filter, s, [4 2], [], stratifyflag, false, latoi, conditions);
    
    statCvsIC  = makesourcecontrast(freq, filter, s, [1 3], [4 2], stratifyflag, true, latoi, conditions);
    statCvsIC2 = makesourcecontrast(freq, filter, s, [1 3], [4 2], stratifyflag, false, latoi, conditions);
    
    % same hemifield contrast congruent minus incongruent
    stat12 = makesourcecontrast(freq, filter, s, [1 2], [], stratifyflag, false, latoi, conditions);
    stat43 = makesourcecontrast(freq, filter, s, [4 3], [], stratifyflag, false, latoi, conditions);
    
    % vs neutral condition (don't stratify RT)
    stat15 = makesourcecontrast(freq, filter, s, [1 5], [], false, false, latoi, conditions);
    stat25 = makesourcecontrast(freq, filter, s, [2 5], [], false, false, latoi, conditions);
    stat35 = makesourcecontrast(freq, filter, s, [3 5], [], false, false, latoi, conditions);
    stat45 = makesourcecontrast(freq, filter, s, [4 5], [], false, false, latoi, conditions);
  else % not yet working FIXME
    % left response
    stat1_12 = makesourcecontrast(freq, filter, s, [1 1 2], [], stratifyflag, false, conditions); % current C, previous C vs IC
    stat1_13 = makesourcecontrast(freq, filter, s, [1 1 3], [], stratifyflag, false, conditions); % current C, previous C vs N
    stat3_45 = makesourcecontrast(freq, filter, s, [3 4 5], [], stratifyflag, false, conditions); % current IC, previous IC vs C
    stat3_56 = makesourcecontrast(freq, filter, s, [3 5 6], [], stratifyflag, false, conditions); % current IC, previous IC vs N
    % right response
    stat4_12 = makesourcecontrast(freq, filter, s, [4 1 2], [], stratifyflag, false, conditions); % current C, previous C vs IC
    stat4_13 = makesourcecontrast(freq, filter, s, [4 1 3], [], stratifyflag, false, conditions); % current C, previous C vs N
    stat2_45 = makesourcecontrast(freq, filter, s, [2 4 5], [], stratifyflag, false, conditions); % current IC, previous IC vs C
    stat2_56 = makesourcecontrast(freq, filter, s, [2 5 6], [], stratifyflag, false, conditions); % current IC, previous IC vs N
  end
end
stat=[];
try stat.stat13=stat13;end
try stat.stat42=stat42;end
try stat.stat12=stat12;end
try stat.stat43=stat43;end
try stat.stat15=stat15;end
try stat.stat25=stat25;end
try stat.stat35=stat35;end
try stat.stat45=stat45;end
try stat.statCvsIC=statCvsIC;end
try stat.statCvsIC2=statCvsIC2;end
try stat.stat1_12=stat1_12;end
try stat.stat1_13=stat1_13;end
try stat.stat3_45=stat3_45;end
try stat.stat3_56=stat3_56;end
try stat.stat4_12=stat4_12;end
try stat.stat4_13=stat4_13;end
try stat.stat2_45=stat2_45;end
try stat.stat2_56=stat2_56;end



function stat = makesourcecontrast(freq, filter, s, contrast, whichflip, stratifyflag, poolhemi, latoi, conditions)

% compute contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest';

if strcmp(latoi, 'post') || (strcmp(latoi, 'pre') && strcmp(conditions, 'previous'))
  tmpcfg            = [];
  if strcmp(latoi, 'post')
    tmpcfg.trials   = find(ismember(freq.trialinfo(:,1),[contrast whichflip]) & freq.trialinfo(:,end)==2); % for the pst trials only
    idx = 1;
  elseif strcmp(latoi, 'pre')
    tmpcfg.trials   = find(ismember(freq.trialinfo(:,end),[contrast whichflip]));
    idx = 4;
  end
  tmpfreq           = ft_selectdata(tmpcfg, freq);
  s.pow             = zeros(numel(s.inside), numel(tmpcfg.trials));
  s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
  s.trialinfo       = tmpfreq.trialinfo;
  
  if ~isempty(whichflip)
    s = fliphemitrials(s, 'pow', idx, whichflip, contrast);
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
  
  cfgs.design                 = nan(1,size(s.trialinfo,1));
  cfgs.design(s.trialinfo(:,idx)==contrast(1)) = 1;
  cfgs.design(s.trialinfo(:,idx)==contrast(2)) = 2;
  stat                        = ft_sourcestatistics(cfgs, s);
  stat                        = rmfield(stat, {'prob', 'cirange', 'mask'});
  
  % save some memory on disk
  if isfield(stat, 'tri'), stat.tri = int16(stat.tri); end
  stat.pos = single(stat.pos);
  
elseif strcmp(latoi, 'pre') && strcmp(conditions, 'current_previous')
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


