function [coh,zx13,zx42,looptime] = vismot_bf_pre_coh(subject,varargin)

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

[freq, tlckpre] =  vismot_spectral_pre(subject,'output','csd','conditions','previous', 'foilim', [frequency frequency], 'smoothing', smoothing);


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
  for k = 1:numel(freq)
    freq(k).grad = ft_transform_geometry(M, ft_convert_units(freq(k).grad, 'm'));
  end
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
allfreq = freq(1);
allfreq.crsspctrm = mean(cat(3,freq.crsspctrm),3);

lambda = 0.1.*trace(allfreq.crsspctrm)./numel(allfreq.label);

if sum(leadfield.inside)>6000
  memreq = 'low';
else
  memreq = 'high';
end

allcoh = estimate_coh2x2_2dip_new(leadfield,allfreq,'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0]);

leadfield_scalar = leadfield;
cnt = 0;
for k = find(leadfield.inside')
  cnt = cnt+1;
  leadfield_scalar.leadfield{k} = leadfield.leadfield{k}*allcoh.ori(cnt,:)';
end
leadfield_scalar = rmfield(leadfield_scalar, 'v');

% this line was to check equivalence of the results.
%cohtmp = estimate_coh2x2_2dip_new(leadfield_scalar,allfreq,'memory','low','lambda',0.05);
for k = 1:5
  coh(k) = estimate_coh2x2_2dip_new(leadfield_scalar,freq(k),'memory',memreq,'lambda',lambda, 'outputflags', [1 0 0 0]);
end

N     = 150;

dcoh13 = single(zeros(size(coh(1).coh)));
dcoh42 = dcoh13;
dcoh13sq = dcoh13;
dcoh42sq = dcoh42;

for k = 1:nrand
  tic;
  indx = sort(randperm(numel(allfreq.label),N)); % keep it sorted!! -> subsampling of sensors
  tmpleadfield = leadfield;
  tmpleadfield.leadfield(tmpleadfield.inside) = cellrowselect(tmpleadfield.leadfield(tmpleadfield.inside),indx);
  tmpcfg.channel = allfreq.label(indx);
  tmpfreq = ft_selectdata(tmpcfg, allfreq);
  tmplambda = 0.1.*trace(tmpfreq.crsspctrm)./numel(tmpfreq.label);
  tmpcoh  = estimate_coh2x2_2dip_new(tmpleadfield,tmpfreq,'memory',memreq,'lambda',tmplambda, 'outputflags', [1 0 0 0]);

  tmpleadfield_scalar = tmpleadfield;
  cnt = 0;
  for m = find(tmpleadfield.inside')
    cnt = cnt+1;
    tmpleadfield_scalar.leadfield{m} = tmpleadfield.leadfield{m}*tmpcoh.ori(cnt,:)';
  end
  tmpleadfield_scalar = rmfield(tmpleadfield_scalar, 'v');
  
  for m = 1:5
    tmpfreq = ft_selectdata(tmpcfg, freq(m));
    tmpcoh(m)  = estimate_coh2x2_2dip_new(tmpleadfield_scalar,tmpfreq,'memory',memreq,'lambda',tmplambda, 'outputflags', [1 0 0 0]);
  end
  
  dcoh13   = dcoh13   + double(abs(tmpcoh(1).coh) - abs(tmpcoh(3).coh));
  dcoh13sq = dcoh13sq + double(abs(tmpcoh(1).coh) - abs(tmpcoh(3).coh)).^2;
  dcoh42   = dcoh42   + double(abs(tmpcoh(4).coh) - abs(tmpcoh(2).coh));
  dcoh42sq = dcoh42sq + double(abs(tmpcoh(4).coh) - abs(tmpcoh(2).coh)).^2;
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
  zx13 = abs(coh(1).coh)-abs(coh(3).coh);
  zx42 = abs(coh(4).coh)-abs(coh(2).coh);
end

if ~exist('looptime', 'var'), looptime = nan; end



%%condition 1: cue left, response left
%%condition 2: cue left, response right
%%condition 3: cue right, response left
%%condition 4: cue right, response right
%%condition 5: catch trial
