function [coh13,coh42,zx13,zx42,looptime] = vismot_bf_coh_stratified(subject,varargin)

%function [source, filter, freq] = vismot_bf_post(subject,varargin)
toi            = ft_getopt(varargin, 'toi', 'pre');
conditions     = ft_getopt(varargin, 'conditions', 'previous');
frequency      = ft_getopt(varargin, 'frequency', 20);
smoothing      = ft_getopt(varargin, 'smoothing', []);
sourcemodel    = ft_getopt(varargin, 'sourcemodel');
nrand          = ft_getopt(varargin, 'nrand', 100); % number of randomization for sensor subsampling
refindx        = ft_getopt(varargin, 'refindx', []);
stratflag      = ft_getopt(varargin, 'stratflag', true);
N              = ft_getopt(varargin, 'N', 75);
lambda         = ft_getopt(varargin, 'lambda', '10%');

if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 12;
  end
end

[freq, tlck] =  vismot_spectral(subject,'toi', toi, 'conditions', conditions, 'output','fourier', 'foilim', [frequency frequency], 'smoothing', smoothing);


% load in the head model and the source model.
if isempty(sourcemodel), error('a sourcemodel should be provided'); end
load(fullfile(subject.pathname, 'headmodel', [subject.name, '_headmodel']));

% coregister the gradiometers if needed
if ~isempty(strfind(subject.datafile, 'h'))
  load(fullfile(subject.pathname, 'dewar2head_avg', [subject.name, 'dewar2head_avg']));
  
  % the transformation matrix M is in centimeters
  M(1:3,4)  = M(1:3,4)./100; % convert to meters
  for k = 1:numel(freq)
    freq(k).grad = ft_transform_geometry(M, ft_convert_units(freq(k).grad, 'm'));
  end
end

headmodel   = ft_convert_units(headmodel,   'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');

% compute leadfields
cfg           = [];
cfg.grad      = freq(1).grad;
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.channel   = freq.label;
cfg.singleshell.batchsize = 2000;
leadfield     = ft_prepare_leadfield(cfg);
for k = find(leadfield.inside')
  tmp = leadfield.leadfield{k};
  [u,s,v] = svd(tmp,'econ');
  leadfield.leadfield{k} = tmp*v(:,1:2);
  leadfield.v{k} = v(:,1:2);
end

if sum(leadfield.inside)>6000
  memreq = 'low';
else
  memreq = 'high';
end

coh13 = estimate_coh2x2_2dip_stratify(leadfield,{freq(1) freq(3)},'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0], 'refindx', refindx, 'stratflag', stratflag);
coh42 = estimate_coh2x2_2dip_stratify(leadfield,{freq(4) freq(2)},'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0], 'refindx', refindx, 'stratflag', stratflag);



dcoh13 = single(zeros(size(coh13.coh)));
dcoh42 = dcoh13;
dcoh13sq = dcoh13;
dcoh42sq = dcoh42;

for k = 1:nrand
  tic;
  indx = sort(randperm(numel(freq(1).label),N)); % keep it sorted!! -> subsampling of sensors
  tmpleadfield = leadfield;
  tmpleadfield.leadfield(tmpleadfield.inside) = cellrowselect(tmpleadfield.leadfield(tmpleadfield.inside),indx);
  
  tmpcfg = [];
  tmpcfg.channel = freq(1).label(indx);
  for m = 1:numel(freq)
    tmpfreq(m) = ft_selectdata(tmpcfg, freq(m));
  end
  
  tmpcoh13 = estimate_coh2x2_2dip_stratify(tmpleadfield,{tmpfreq(1) tmpfreq(3)},'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0], 'refindx', refindx, 'stratflag', stratflag);
  tmpcoh42 = estimate_coh2x2_2dip_stratify(tmpleadfield,{tmpfreq(4) tmpfreq(2)},'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0], 'refindx', refindx, 'stratflag', stratflag);


  dcoh13   = dcoh13   + double(abs(tmpcoh13.coh_1) - abs(tmpcoh13.coh_2));
  dcoh13sq = dcoh13sq + double(abs(tmpcoh13.coh_1) - abs(tmpcoh13.coh_2)).^2;
  dcoh42   = dcoh42   + double(abs(tmpcoh42.coh_1) - abs(tmpcoh42.coh_2));
  dcoh42sq = dcoh42sq + double(abs(tmpcoh42.coh_1) - abs(tmpcoh42.coh_2)).^2;
  looptime(k) = toc;
end

if nrand>0
  mx13 = dcoh13./nrand;
  sx13 = sqrt((dcoh13sq-(dcoh13.^2)./nrand)./nrand);
  zx13 = mx13./sx13;
  
  mx42 = dcoh42./nrand;
  sx42 = sqrt((dcoh42sq-(dcoh42.^2)./nrand)./nrand);
  zx42 = mx42./sx42;
else
  zx13 = abs(coh13.coh_1)-abs(coh13.coh_2);
  zx42 = abs(coh42.coh_1)-abs(coh42.coh_2);
end

if ~exist('looptime', 'var'), looptime = nan; end



%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial
