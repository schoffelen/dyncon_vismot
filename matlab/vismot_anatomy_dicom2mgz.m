function vismot_anatomy_dicom2mgz

subjinfo;

for k = 5:numel(SUBJ)
  if exist(fullfile(SUBJ(k).pathname,'mri',[SUBJ(k).name,'.mgz']),'file')
    continue;
  end
  
  if ~isempty(SUBJ(k).mriname)
    name = SUBJ(k).mriname;
  else
    name = SUBJ(k).name;
  end
  fprintf('processing subject %s\n',SUBJ(k).name);

  cd(SUBJ(k).rawpath);
  cd([name,'_RAW']);
  d   = dir([upper(name),'*']);
  mri = ft_read_mri(d(end).name);
  
  cfg = [];
  cfg.filename = fullfile(SUBJ(k).pathname,'mri',[SUBJ(k).name,'.mgz']);
  cfg.filetype = 'mgz';
  cfg.parameter = 'anatomy';
  ft_volumewrite(cfg, mri);
  clear mri;
  
end
