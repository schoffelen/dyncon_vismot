function [source, parcellation, filter_avgpos, filter_lfsvd] = vismot_bf_lcmv_pre(subject, tlckpre, varargin)

truncate   = ft_getopt(varargin, 'truncate', []);

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
[leadfield, leadfieldorig] = vismot_forward_parcellate(subject, tlckpre);
	
for k = 1:size(leadfield.pos,1)
  if ~isempty(truncate)
	   tmps = leadfield.s{k};
		 tmps = cumsum(tmps)./sum(tmps);
		 leadfield.leadfield{k} = leadfield.u{k}(:,tmps<truncate);
	else
		leadfield.leadfield{k,1} = leadfield.u{k};
	end
end

tlck = tlckpre(1);
tlck.cov = (tlckpre(1).cov+tlckpre(2).cov+tlckpre(3).cov+tlckpre(4).cov+tlckpre(5).cov)./5;

cfg                 = [];
cfg.grid            = rmfield(leadfield, 'leadfield');
cfg.headmodel       = headmodel;
cfg.method          = 'lcmv';
cfg.keeptrials      = 'yes';
cfg.lcmv.lambda     = '100%';
cfg.lcmv.keepfilter = 'yes';
source1             = ft_sourceanalysis(cfg, tlck);
filter_avgpos       = source1.avg.filter;
cfg.grid            = leadfield;
source2             = ft_sourceanalysis(cfg, tlck);
filter_lfsvd        = source2.avg.filter;
cfg.grid            = leadfieldorig;
source              = ft_sourceanalysis(cfg, tlck);

load(fullfile('/project/3011085.03/analysis/mri/Conte69_32k/atlas_subparc374_8k'));
[source, parcellation] = vismot_bf_lcmv_parcellate(source, tlck, 'parcellation', atlas);
