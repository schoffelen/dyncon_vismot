function headmodel = vismot_anatomy_headmodel(subject)

if nargin<1
  % do it for all subjects
  subjects = vismot_subjinfo;
  for k = [1:numel(subjects)]
    headmodel = vismot_anatomy_headmodel(subjects(k).name);
  end
  return;
end

if ischar(subject)
  subject = vismot_subjinfo(subject);
end

pname = fullfile(subject.pathname,'mri');
mri   = ft_read_mri(fullfile(pname,[subject.name,'.mgz']));
load(fullfile(pname,[subject.name,'_transform_vox2bti.mat']));
mri.transform = transform;
mri.coordsys  = 'bti';

cfg = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 10000;
bnd = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, bnd);

pname = fullfile(subject.pathname,'headmodel');
fname = fullfile(pname, [subject.name,'_headmodel']);
save(fname, 'headmodel');


