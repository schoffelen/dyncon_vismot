function vismot_anatomy_skullstrip

subjinfo;

for k = [1:3 5:numel(SUBJ)]
  fname = fullfile(SUBJ(k).pathname,'mri',SUBJ(k).name);
  mri   = ft_read_mri([fname,'.mgz']);

  threshold = 0.5;
  T         = inv(mri.transform);
  center    = round(T(1:3,4))';
  subjectname = SUBJ(k).name;

  % save temporarily as a nii
  t   = tempname;
  
  cfg = [];
  cfg.filename = t;
  cfg.filetype = 'nifti';
  cfg.parameter = 'anatomy';
  ft_volumewrite(cfg, mri);
  
  str = ['/opt/fsl/5.0.9/bin/bet ',t,'.nii ',t];
  str = [str,'-R -f ',num2str(threshold),' -c ', num2str(center),' -g 0 -m -v'];
  
  system(str);
  seg  = ft_read_mri([t,'-R.nii.gz']);
  delete([t,'.nii']);
  delete([t,'-R.nii.gz']);
  delete([t,'-R_mask.nii.gz']);
  
  cfg = [];
  cfg.filename = [fname,'skullstrip'];
  cfg.filetype = 'mgz';
  cfg.parameter = 'anatomy';
  ft_volumewrite(cfg, seg);
end
