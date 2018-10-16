function [source, stat13, stat42] = vismot_bf_post(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency = ft_getopt(varargin, 'frequency', 10);
smoothing = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 12;
  end
end

[freqpre,freqpst] = vismot_spectral_prepost(subject,'foilim',[frequency frequency],'smoothing',smoothing,'output','fourier');

for k = 1:5
  % add marker for pre/pst
  freqpre(k).trialinfo(:,end+1) = 1;
  freqpst(k).trialinfo(:,end+1) = 2;
end

cfg = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq = ft_appendfreq(cfg, freqpre(1), freqpre(2), freqpre(3), freqpre(4), freqpre(5), freqpst(1), freqpst(2), freqpst(3), freqpst(4), freqpst(5));
clear freqpst;

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
end

%sourcemodel.inside(11:end)=false;

% compute beamformer common spatial filters
cfg           = [];
cfg.grad      = freq.grad;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.channel   = freq.label;
cfg.singleshell.batchsize = 2000;
leadfield     = ft_prepare_leadfield(cfg);

cfg                 = [];
cfg.grid            = leadfield;
cfg.headmodel       = headmodel;
cfg.method          = 'dics';
cfg.keeptrials      = 'yes';
cfg.dics.lambda     = '10%';
cfg.dics.fixedori   = 'yes';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfilter = 'yes';
tmpsource = ft_sourceanalysis(cfg, freq);
filter    = tmpsource.avg.filter;

cfg2             = [];
cfg2.fwhm        = 'yes';
if ~isfield(sourcemodel, 'dim'), cfg2.fwhmmethod  = 'gaussfit'; end
cfg2.fwhmmaxdist = 0.02;
fwhm             = ft_sourcedescriptives(cfg2, tmpsource);
fwhm             = fwhm.fwhm;

cfg.grid.filter     = filter;
cfg.dics.keepfilter = 'no';

% compute 1-3 and 2-4 contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest';

s     = keepfields(tmpsource,{'freq' 'tri' 'inside' 'pos'});

cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,1),[1 3]) & freq.trialinfo(:,end)==2); % for the pst trials only
tmpfreq           = ft_selectdata(cfg2, freq);
s.pow = zeros(numel(s.inside),numel(cfg2.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
s.trialinfo       = tmpfreq.trialinfo;

cfgs.design                 = s.trialinfo(:,1)';
cfgs.design(cfgs.design==3) = 2;
stat13                      = ft_sourcestatistics(cfgs, s);
stat13 = rmfield(stat13, {'prob', 'cirange', 'mask'});
try, stat13.tri = int16(stat13.tri); end
stat13.pos = single(stat13.pos);

cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,1),[2 4]) & freq.trialinfo(:,end)==2); % for the pst trials only
tmpfreq           = ft_selectdata(cfg2, freq);
s.pow = zeros(numel(s.inside),numel(cfg2.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
s.trialinfo       = tmpfreq.trialinfo;

cfgs.design                 = s.trialinfo(:,1)';
cfgs.design(cfgs.design==4) = 1;
stat42                      = ft_sourcestatistics(cfgs, s);
stat42 = rmfield(stat42, {'prob', 'cirange', 'mask'});
try, stat42.tri = int16(stat42.tri); end
stat42.pos = single(stat42.pos);

try
% compute condition specific power
for k = 1:5
  cfg2.trials = find(freq.trialinfo(:,1)==k & freq.trialinfo(:,end)==2); % for the pst trials only
  tmp         = ft_sourceanalysis(cfg, ft_selectdata(cfg2, freq));
  tmp.fwhm    = fwhm;
  tmp         = smooth_source(tmp, 'parameter', 'pow', 'maxdist', 0.025);
  source(k)   = tmp;
end
end

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial