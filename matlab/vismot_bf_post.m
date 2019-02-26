function [source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC, zx13, zx42] = vismot_bf_post(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency = ft_getopt(varargin, 'frequency', 10);
smoothing = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');
prewhiten = istrue(ft_getopt(varargin, 'prewhiten', false));
lambda = ft_getopt(varargin, 'lambda', ' 100%');
nrand          = ft_getopt(varargin, 'nrand', 0); % number of randomization for sensor subsampling
N              = ft_getopt(varargin, 'N', 90);


if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
  end
end

%this is the old version of the spectral function
%[freqpre,freqpst] = vismot_spectral_prepost(subject,'foilim',[frequency frequency],'smoothing',smoothing,'output','fourier');

freq    = vismot_spectral(subject, 'foilim', [frequency frequency], 'toi', 'prepost', 'balance', true, 'smoothing', smoothing, 'output', 'fourier', 'prewhiten', prewhiten);
freqpre = freq(1,:);
freqpst = freq(2,:);
clear freq;

for k = 1:5
  % add marker for pre/pst
  freqpre(k).trialinfo(:,end+1) = 1;
  freqpst(k).trialinfo(:,end+1) = 2;
end

cfg = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq = ft_appendfreq(cfg, freqpre(1), freqpre(2), freqpre(3), freqpre(4), freqpre(5), freqpst(1), freqpst(2), freqpst(3), freqpst(4), freqpst(5));

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

[source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = sourcepow_post(freq, headmodel, leadfield, lambda);

dpow13 = single(zeros(size(stat13.stat)));
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
  
  [~, tmppow13, tmppow42] = sourcepow_post(tmpfreq, headmodel, tmpleadfield, lambda, true);
  
  dpow13   = dpow13   + tmppow13.stat;
  dpow13sq = dpow13sq + tmppow13.stat.^2;
  dpow42   = dpow42   + tmppow42.stat;
  dpow42sq = dpow42sq + tmppow42.stat.^2;
  looptime(k) = toc;
end

if nrand>0
  mx13 = dpow13./nrand;
  sx13 = sqrt((dpow13sq-(dpow13.^2)./nrand)./nrand);
  zx13 = mx13./sx13;
  
  mx42 = dpow42./nrand;
  sx42 = sqrt((dpow42sq-(dpow42.^2)./nrand)./nrand);
  zx42 = mx42./sx42;
else
  zx13 = [];
  zx42 = [];
end


function [source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = sourcepow_post(freq, headmodel, leadfield, lambda, onlycompute_stat13_42)
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

s     = keepfields(tmpsource, {'freq' 'tri' 'inside' 'pos' 'dim'});


% same response hand contrast congruent minus incongruent
stat13 = makesourcecontrast(freq, filter, s, [1 3], [], false, false);
stat42 = makesourcecontrast(freq, filter, s, [4 2], [], false, false);

if ~onlycompute_stat13_42
statCvsIC  = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, true);
statCvsIC2 = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, false);

% same hemifield contrast congruent minus incongruent
stat12 = makesourcecontrast(freq, filter, s, [1 2], [], false, false);
stat43 = makesourcecontrast(freq, filter, s, [4 3], [], false, false);

% vs neutral condition (don't stratify RT)
stat15 = makesourcecontrast(freq, filter, s, [1 5], [], false, false);
stat25 = makesourcecontrast(freq, filter, s, [2 5], [], false, false);
stat35 = makesourcecontrast(freq, filter, s, [3 5], [], false, false);
stat45 = makesourcecontrast(freq, filter, s, [4 5], [], false, false);

% compute condition specific power, this is without stratification for RT
cfg.grid.filter     = filter;
cfg.dics.keepfilter = 'no';
for k = 1:5
  tmpcfg.trials = find(freq.trialinfo(:,1)==k & freq.trialinfo(:,end)==2); % for the pst trials only
  source(k)     = ft_sourceanalysis(cfg, ft_selectdata(tmpcfg, freq));
end
else
  statCvsIC=[];
  stat15 = [];
  stat25 = [];
  stat35 = [];
  stat45 = [];
  source = [];
end



function stat = makesourcecontrast(freq, filter, s, contrast, whichflip, stratifyflag, poolhemi)

% compute contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest'; 
%cfgs.statistic       = 'indepsamplesT';

tmpcfg            = [];
tmpcfg.trials     = find(ismember(freq.trialinfo(:,1),[contrast whichflip]) & freq.trialinfo(:,end)==2); % for the pst trials only
tmpfreq           = ft_selectdata(tmpcfg, freq);
s.pow             = zeros(numel(s.inside), numel(tmpcfg.trials));
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

cfgs.design                 = nan(1,size(s.trialinfo,1));
cfgs.design(s.trialinfo(:,1)==contrast(1)) = 1;
cfgs.design(s.trialinfo(:,1)==contrast(2)) = 2;
stat                        = ft_sourcestatistics(cfgs, s);
stat                        = rmfield(stat, {'prob', 'cirange', 'mask'});

% save some memory on disk
if isfield(stat, 'tri'), stat.tri = int16(stat.tri); end
stat.pos = single(stat.pos);
