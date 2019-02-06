function [source, stat] = vismot_bf_pre(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)

frequency   = ft_getopt(varargin, 'frequency', 20);
smoothing   = ft_getopt(varargin, 'smoothing', []);
sourcemodel = ft_getopt(varargin, 'sourcemodel');
conditions  = ft_getopt(varargin, 'conditions', 'previous');
prewhiten   = ft_getopt(varargin, 'prewhiten', false);
if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
  end
end

[freq, tlckpre] =  vismot_spectral(subject,'output','fourier','conditions',conditions, 'foilim', [frequency frequency], 'smoothing', smoothing, 'prewhiten', prewhiten);
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



% compute condition specific power
for k = 1:nconditions
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



% compute contrasts as a yuen-welch T value
cd /project/3011085.03/scripts/dyncon_vismot/matlab/private/
cfgs        = [];
cfgs.method = 'montecarlo';
cfgs.numrandomization = 0;
cfgs.statistic        = 'statfun_yuenTtest'; % This statistics function
% is not available.
%cfgs.statistic       = 'indepsamplesT';

s     = keepfields(tmpsource,{'freq' 'tri' 'inside' 'pos' 'dim'});

switch conditions
  case 'previous'
    % same response hand contrast
    % left hand response, C vs neutral
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[1 5]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==5) = 2;
    stat15    = ft_sourcestatistics(cfgs, s);
    stat15 = rmfield(stat15, {'prob', 'cirange', 'mask'});
    try, stat15.tri = int16(stat15.tri); end
    stat15.pos = single(stat15.pos);
    
    % left hand response, IC vs neutral
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[3 5]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==3) = 1;
    cfgs.design(cfgs.design==5) = 2;
    stat35    = ft_sourcestatistics(cfgs, s);
    stat35 = rmfield(stat35, {'prob', 'cirange', 'mask'});
    try, stat35.tri = int16(stat35.tri); end
    stat35.pos = single(stat35.pos);
    
    % right hand response, C vs neutral
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[4 5]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==4) = 1;
    cfgs.design(cfgs.design==5) = 2;
    stat45    = ft_sourcestatistics(cfgs, s);
    stat45 = rmfield(stat45, {'prob', 'cirange', 'mask'});
    try, stat45.tri = int16(stat45.tri); end
    stat45.pos = single(stat45.pos);
    
    % right hand response, IC vs neutral
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[2 5]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==2) = 1;
    cfgs.design(cfgs.design==5) = 2;
    stat25    = ft_sourcestatistics(cfgs, s);
    stat25 = rmfield(stat25, {'prob', 'cirange', 'mask'});
    try, stat25.tri = int16(stat25.tri); end
    stat25.pos = single(stat25.pos);
    
    % left hand response, C-IC
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[1 3]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==3) = 2;
    stat13                      = ft_sourcestatistics(cfgs, s);
    stat13 = rmfield(stat13, {'prob', 'cirange', 'mask'});
    try, stat13.tri = int16(stat13.tri); end
    stat13.pos = single(stat13.pos); % what is this step for?
    
    % right hand response, C-IC
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,end),[2 4]));
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,4)';
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
    
    cfgs.design                 = s.trialinfo(:,4)';
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
    
    cfgs.design                 = s.trialinfo(:,4)';
    cfgs.design(cfgs.design==4) = 1;
    cfgs.design(cfgs.design==3) = 2;
    stat43                      = ft_sourcestatistics(cfgs, s);
    stat43 = rmfield(stat43, {'prob', 'cirange', 'mask'});
    try, stat43.tri = int16(stat43.tri); end
    stat43.pos = single(stat43.pos);
    
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
    
    % very ugly. FIXME
    stat.stat13 = stat13;
    stat.stat42 = stat42;
    stat.stat12 = stat12;
    stat.stat43 = stat43;
    stat.stat15 = stat15;
    stat.stat25 = stat25;
    stat.stat35 = stat35;
    stat.stat45 = stat45;
    
  case 'current_previous'
    % left response
    % current C, previous C vs IC
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[1 2]) & freq.trialinfo(:,1)==1);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    stat1    = ft_sourcestatistics(cfgs, s);
    stat1 = rmfield(stat1, {'prob', 'cirange', 'mask'});
    try, stat1.tri = int16(stat1.tri); end
    stat1.pos = single(stat1.pos);
    
    % current C, previous C vs N
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[1 3]) & freq.trialinfo(:,1)==1);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==3) = 2;
    stat2    = ft_sourcestatistics(cfgs, s);
    stat2 = rmfield(stat2, {'prob', 'cirange', 'mask'});
    try, stat2.tri = int16(stat2.tri); end
    stat2.pos = single(stat2.pos);
    
    % current IC, previous IC vs C
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[4 5]) & freq.trialinfo(:,1)==3);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==5) = 1;
    cfgs.design(cfgs.design==4) = 2;
    stat3    = ft_sourcestatistics(cfgs, s);
    stat3 = rmfield(stat3, {'prob', 'cirange', 'mask'});
    try, stat3.tri = int16(stat3.tri); end
    stat3.pos = single(stat3.pos);
    
    % current IC, previous IC vs N
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[5 6]) & freq.trialinfo(:,1)==3);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==5) = 1;
    cfgs.design(cfgs.design==6) = 2;
    stat4    = ft_sourcestatistics(cfgs, s);
    stat4 = rmfield(stat4, {'prob', 'cirange', 'mask'});
    try, stat4.tri = int16(stat4.tri); end
    stat4.pos = single(stat4.pos);
    
    
    % right response
    % current C, previous C vs IC
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[1 2]) & freq.trialinfo(:,1)==4);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    stat5    = ft_sourcestatistics(cfgs, s);
    stat5 = rmfield(stat5, {'prob', 'cirange', 'mask'});
    try, stat5.tri = int16(stat5.tri); end
    stat5.pos = single(stat5.pos);
    
    % current C, previous C vs N
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[1 3]) & freq.trialinfo(:,1)==4);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==3) = 2;
    stat6    = ft_sourcestatistics(cfgs, s);
    stat6 = rmfield(stat6, {'prob', 'cirange', 'mask'});
    try, stat6.tri = int16(stat6.tri); end
    stat6.pos = single(stat6.pos);
    
    % current IC, previous IC vs C
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[4 5]) & freq.trialinfo(:,1)==2);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==5) = 1;
    cfgs.design(cfgs.design==4) = 2;
    stat7    = ft_sourcestatistics(cfgs, s);
    stat7 = rmfield(stat7, {'prob', 'cirange', 'mask'});
    try, stat7.tri = int16(stat7.tri); end
    stat7.pos = single(stat7.pos);
    
    % current IC, previous IC vs N
    cfg2              = [];
    cfg2.trials       = find(ismember(freq.trialinfo(:,5),[5 6]) & freq.trialinfo(:,1)==2);
    tmpfreq           = ft_selectdata(cfg2, freq);
    s.pow = zeros(numel(s.inside),numel(cfg2.trials));
    s.pow(s.inside,:) = fourier2pow(cat(3, filter{:}), tmpfreq.fourierspctrm, tmpfreq.cumtapcnt);
    s.trialinfo       = tmpfreq.trialinfo;
    
    cfgs.design                 = s.trialinfo(:,end)';
    cfgs.design(cfgs.design==5) = 1;
    cfgs.design(cfgs.design==6) = 2;
    stat8    = ft_sourcestatistics(cfgs, s);
    stat8 = rmfield(stat8, {'prob', 'cirange', 'mask'});
    try, stat8.tri = int16(stat8.tri); end
    stat8.pos = single(stat8.pos);
    
    stat.stat1=stat1;
    stat.stat2=stat2;
    stat.stat3=stat3;
    stat.stat4=stat4;
    stat.stat5=stat5;
    stat.stat6=stat6;
    stat.stat7=stat7;
    stat.stat8=stat8;
end




%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial