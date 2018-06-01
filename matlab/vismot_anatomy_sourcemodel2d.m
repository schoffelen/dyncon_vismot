function sourcemodel = vismot_anatomy_sourcemodel2d(subject)

% VISMOT_ANATOMY_SOURCEMODEL2D loads in the workbench files (8k) and
% coregisters the cortical sheet to the MEG coordinate system

% load in the cortical sheet
datapath = fullfile(subject.pathname,'mri',subject.name,'workbench');
filename = fullfile(datapath,[subject.name,'.L.midthickness.8k_fs_LR.surf.gii']);
sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

% get the necessary coregistration information
datapath = fullfile(subject.pathname,'mri');
load(fullfile(datapath,[subject.name,'_transform_vox2mni']));
T1 = transform;
load(fullfile(datapath,[subject.name,'_transform_vox2bti']));
T2 = transform;

sourcemodel = ft_transform_geometry((T2/T1), sourcemodel);
sourcemodel.inside = sourcemodel.atlasroi>0;
sourcemodel = rmfield(sourcemodel, 'atlasroi');

