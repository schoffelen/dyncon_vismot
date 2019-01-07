function [source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = vismot_bf_post(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency = ft_getopt(varargin, 'frequency', 10);
smoothing = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
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

tmpcfg             = [];
tmpcfg.fwhm        = 'yes';
if ~isfield(sourcemodel, 'dim')
  tmpcfg.fwhmmethod  = 'gaussfit';
  tmpcfg.fwhmmaxdist = 0.02;
end
fwhm             = ft_sourcedescriptives(tmpcfg, tmpsource);
fwhm             = fwhm.fwhm;

cfg.grid.filter     = filter;
cfg.dics.keepfilter = 'no';


s     = keepfields(tmpsource, {'freq' 'tri' 'inside' 'pos' 'dim'});

statCvsIC = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, true);
statCvsIC2 = makesourcecontrast(freq, filter, s, [1 3], [4 2], false, false);


% same response hand contrast congruent minus incongruent
stat13 = makesourcecontrast(freq, filter, s, [1 3], [], false, false);
stat42 = makesourcecontrast(freq, filter, s, [4 2], [], false, false);

% same hemifield contrast congruent minus incongruent
stat12 = makesourcecontrast(freq, filter, s, [1 2], [], false, false);
stat43 = makesourcecontrast(freq, filter, s, [4 3], [], false, false);

% vs neutral condition (don't stratify RT)
stat15 = makesourcecontrast(freq, filter, s, [1 5], [], false, false);
stat25 = makesourcecontrast(freq, filter, s, [2 5], [], false, false);
stat35 = makesourcecontrast(freq, filter, s, [3 5], [], false, false);
stat45 = makesourcecontrast(freq, filter, s, [4 5], [], false, false);

% compute condition specific power, this is without stratification for RT
for k = 1:5
  tmpcfg.trials = find(freq.trialinfo(:,1)==k & freq.trialinfo(:,end)==2); % for the pst trials only
  tmp         = ft_sourceanalysis(cfg, ft_selectdata(tmpcfg, freq));
%   %try
%     tmp.fwhm    = fwhm;
%     tmp.inside  = tmp.inside & isfinite(fwhm);
%     tmp         = smooth_source(tmp, 'parameter', 'pow', 'maxdist', 0.025);
%   %catch
%   %  tmp = removefields(tmp, 'fwhm');
%   %end
  source(k)   = tmp;
end

% % smooth contrasts
% % same response hand:
% tmp         = stat13;
% tmp.fwhm    = fwhm;
% tmp.inside  = tmp.inside&isfinite(fwhm);
% tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
% stat13      = tmp;
% 
% tmp         = stat42;
% tmp.fwhm    = fwhm;
% tmp.inside  = tmp.inside&isfinite(fwhm);
% tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
% stat42         = tmp;
% 
% % same hemifield
% tmp         = stat12;
% tmp.fwhm    = fwhm;
% tmp.inside  = tmp.inside&isfinite(fwhm);
% tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
% stat12      = tmp;
% 
% tmp         = stat43;
% tmp.fwhm    = fwhm;
% tmp.inside  = tmp.inside&isfinite(fwhm);
% tmp         = smooth_source(tmp, 'parameter', 'stat', 'maxdist', 0.025);
% stat43      = tmp;

%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial

function stat = makesourcecontrast(freq, filter, s, contrast, whichflip, stratifyflag, poolhemi)

% compute contrasts as a yuen-welch T value
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest'; % This statistics function
% is not available.
%cfgs.statistic       = 'indepsamplesT';

tmpcfg            = [];
tmpcfg.trials     = find(ismember(freq.trialinfo(:,1),[contrast whichflip]) & freq.trialinfo(:,end)==2); % for the pst trials only
tmpfreq           = ft_selectdata(tmpcfg, freq);
s.pow             = zeros(numel(s.inside), numel(tmpcfg.trials));
s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);  
s.trialinfo       = tmpfreq.trialinfo;

if ~isempty(whichflip)
    s = fliphemitrials(s, 'pow', 1, [whichflip], [contrast]);
end
if poolhemi
    load standard_sourcemodel3d4mm;
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
