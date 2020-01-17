function vismot_anatomy_dicom2mgz(subject)

if nargin<1
  subject = vismot_subjinfo;
elseif ischar(subject)
  subject = vismot_subjinfo(subject);
end
  

for k = 1:numel(subject)
  if exist(fullfile(subject(k).pathname,'mri',[subject(k).name,'.mgz']),'file')
    continue;
  end
  
  if isfield(subject(k), 'mriname') && ~isempty(subject(k).mriname)
    name = subject(k).mriname;
  else
    name = subject(k).name;
  end
  fprintf('processing subject %s\n',subject(k).name);

  cd(subject(k).rawpath);
  cd([name,'/anat']);
  d   = dir([upper(name),'*']);
  mri = ft_read_mri(d(end).name);
  
  cfg = [];
  cfg.filename = fullfile(subject(k).pathname,'mri',[subject(k).name,'.mgz']);
  cfg.filetype = 'mgz';
  cfg.parameter = 'anatomy';
  ft_volumewrite(cfg, mri);
  clear mri;
  
end
