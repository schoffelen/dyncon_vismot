function [source, stat13, stat42, stat12, stat43] = vismot_bf_pre(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency   = ft_getopt(varargin, 'frequency', 20);
smoothing   = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');
nrand       = ft_getopt(varargin, 'nrand', 100); % number of randomization for sensor subsampling

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 12;
  end
end

[freq, tlckpre] =  vismot_spectral_pre(subject,'output','fourier','conditions','previous', 'foilim', [frequency frequency], 'smoothing', smoothing);
for k = 1:numel(freq)
  if ~isfield(freq(k),'trialinfo')
    freq(k).trialinfo = ones(numel(freq(k).cumtapcnt),1).*k;
  else
    freq(k).trialinfo(:,end+1) = k;
  end
end

cfg = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq = ft_appendfreq(cfg, freq(1), freq(2), freq(3), freq(4), freq(5));
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
cfg.dics.lambda     = '100%';
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

cfg.grid.filter     = filter;
cfg.dics.keepfilter = 'no';

% compute 1-3 and 2-4 contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest'; % This statistics function
% is not available.
%cfgs.statistic       = 'indepsamplesT';

s     = keepfields(tmpsource,{'freq' 'tri' 'inside' 'pos' 'dim'});

% same response hand contrast
cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,end),[1 3]));
tmpfreq           = ft_selectdata(cfg2, freq);
s.pow = zeros(numel(s.inside),numel(cfg2.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);  
s.trialinfo       = tmpfreq.trialinfo;

cfgs.design                 = s.trialinfo(:,1)';
cfgs.design(cfgs.design==3) = 2;
stat13                      = ft_sourcestatistics(cfgs, s);
stat13 = rmfield(stat13, {'prob', 'cirange', 'mask'});
try, stat13.tri = int16(stat13.tri); end
stat13.pos = single(stat13.pos); % what is this step for?

cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,end),[2 4]));
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

% same hemifield contrast
cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,end),[1 2]));
tmpfreq           = ft_selectdata(cfg2, freq);
s.pow = zeros(numel(s.inside),numel(cfg2.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
s.trialinfo       = tmpfreq.trialinfo;

cfgs.design                 = s.trialinfo(:,1)';
stat12                      = ft_sourcestatistics(cfgs, s);
stat12 = rmfield(stat12, {'prob', 'cirange', 'mask'});
try, stat12.tri = int16(stat12.tri); end
stat12.pos = single(stat12.pos);

cfg2              = [];
cfg2.trials       = find(ismember(freq.trialinfo(:,end),[4 3]));
tmpfreq           = ft_selectdata(cfg2, freq);
s.pow = zeros(numel(s.inside),numel(cfg2.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
s.trialinfo       = tmpfreq.trialinfo;

cfgs.design                 = s.trialinfo(:,1)';
cfgs.design(cfgs.design==4) = 1;
cfgs.design(cfgs.design==3) = 2;
stat43                      = ft_sourcestatistics(cfgs, s);
stat43 = rmfield(stat43, {'prob', 'cirange', 'mask'});
try, stat43.tri = int16(stat43.tri); end
stat43.pos = single(stat43.pos);

% compute condition specific power
for k = 1:5
  cfg2.trials = find(freq.trialinfo(:,end)==k);
  tmp         = ft_sourceanalysis(cfg, ft_selectdata(cfg2, freq));
%   %try
%     tmp.fwhm    = fwhm;
%     tmp.inside  = tmp.inside & isfinite(fwhm);
%     tmp         = smooth_source(tmp, 'parameter', 'pow', 'maxdist', 0.025);
%   %catch
%   %  tmp = removefields(tmp, 'fwhm');
%   %end
  source(k)   = tmp;
end

% smooth contrasts
% same response hand:
tmp         = stat13;
tmp.fwhm    = fwhm;
tmp.inside  = tmp.inside&isfinite(fwhm);
%tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
stat13      = tmp;

tmp         = stat42;
tmp.fwhm    = fwhm;
tmp.inside  = tmp.inside&isfinite(fwhm);
%tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
stat42         = tmp;

% same hemifield
tmp         = stat12;
tmp.fwhm    = fwhm;
tmp.inside  = tmp.inside&isfinite(fwhm);
%tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
stat12      = tmp;

tmp         = stat43;
tmp.fwhm    = fwhm;
tmp.inside  = tmp.inside&isfinite(fwhm);
%tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
stat43      = tmp;

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial