function [leadfield, leadfieldorig] = vismot_forward_parcellate(subject, tlckpre, varargin)

if ischar(subject)
	subject = vismot_subjinfo(subject);
end

% load in the head model and the source model.
sourcemodel = vismot_anatomy_sourcemodel2d(subject);
load(fullfile(subject.pathname, 'headmodel', [subject.name, '_headmodel']));

% coregister the gradiometers if needed
if ~isempty(strfind(subject.datafile, 'h'))
  load(fullfile(subject.pathname, 'dewar2head_avg', [subject.name, 'dewar2head_avg']));
  
  % the transformation matrix M is in centimeters
  M(1:3,4) = M(1:3,4)./100; % convert to meters
  grad     = ft_transform_geometry(M, ft_convert_units(tlckpre(1).grad, 'm'));
else
	grad = ft_convert_units(tlckpre(1).grad, 'm');
end

headmodel   = ft_convert_units(headmodel, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');
if ~isfield(sourcemodel, 'inside')
  sourcemodel.inside = true(size(sourcemodel.pnt,1),1);
end

%sourcemodel.inside(11:end)=false;

% compute beamformer common spatial filters
cfg = [];
cfg.grad      = grad;
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.channel = tlckpre(1).label;
cfg.backproject = 'no';
cfg.singleshell.batchsize = 2000;
leadfieldorig = ft_prepare_leadfield(cfg);

% load the atlas
a = load(fullfile('/project/3011085.03/analysis/mri/Conte69_32k/atlas_subparc374_8k')); % load into structure, otherwise might be mistaken by function 'atlas'.

leadfield = rmfield(leadfieldorig, {'pos' ,'tri', 'inside', 'unit', 'leadfield'});
leadfield.brainordinate = a.atlas;
leadfield.brainordinate.pos = leadfieldorig.pos;
leadfield.brainordinate.tri = leadfieldorig.tri;

for k = 1:max(a.atlas.parcellation)
	sel = a.atlas.parcellation==k;
  lf  = cat(2,leadfieldorig.leadfield{sel});
	[u,s,v] = svd(lf,'econ');
	pos = mean(leadfield.brainordinate.pos(sel,:));
	leadfield.pos(k,:) = pos;
	leadfield.u{k,1} = u;
	leadfield.s{k,1} = diag(s);
	leadfield.v{k,1} = v;
end
