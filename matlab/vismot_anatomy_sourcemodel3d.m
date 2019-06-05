function [sourcemodel] = vismot_anatomy_sourcemodel3d(subject, resolution)

% VISMOT_ANATOMY_SOURCEMODEL3D computes a 3D regular grid with 
% specified resolution based on an inverse warp of a template
% grid in MNI space.

% load the mri + coreg
mri = ft_read_mri(fullfile(subject.pathname,'mri',subject.name,'mri',sprintf('%s.mgz',subject.name)));
tmp = load(fullfile(subject.pathname,'mri',sprintf('%s_transform_vox2bti.mat',subject.name)));
mri.transform = tmp.transform;
mri.coordsys  = 'bti';

% create the grid
cfg = [];
cfg.grid.warpmni    = 'yes';
cfg.grid.resolution = resolution;
cfg.grid.nonlinear  = 'yes';
cfg.mri = mri;
sourcemodel = ft_prepare_sourcemodel(cfg);

% remove the mri-structure from grid.cfg
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'mri');
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'callinfo');

save(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d%dmm',subject.name, resolution)),'sourcemodel');

