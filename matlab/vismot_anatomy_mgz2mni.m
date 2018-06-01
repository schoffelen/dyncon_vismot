function vismot_anatomy_mgz2mni

subjinfo;

for k = [1:3 5:numel(SUBJ)]
  if exist(fullfile(SUBJ(k).pathname,'mri',[SUBJ(k).name,'_transform_vox2mni.mat']),'file')
    fname = fullfile(SUBJ(k).pathname,'mri',SUBJ(k).name);
    mri   = ft_read_mri([fname,'.mgz']);
    
    load(fullfile(SUBJ(k).pathname,'mri',[SUBJ(k).name,'_transform_vox2mni.mat']));
    if isequal(mri.transform, transform)
      % do nothing
    else
      fprintf('adding coregistration information to the mrifile of subject %s\n', SUBJ(k).name);
      
      mri.transform = transform;
      
      cfg = [];
      cfg.parameter = 'anatomy';
      cfg.filename  = fname;
      cfg.filetype  = 'mgz';
      ft_volumewrite(cfg, mri);
      
      clear mri;
    end
    
  else
    fname = fullfile(SUBJ(k).pathname,'mri',SUBJ(k).name);
    mri   = ft_read_mri([fname,'.mgz']);
    
    cfg = [];
    cfg.coordsys = 'spm';
    mri = ft_volumerealign(cfg, mri);
    
    transform = mri.transform;
    mriname   = [fname,'.mgz'];
    save([fname,'_transform_vox2mni'],'transform','mriname');
    
    clear mri;
  end  
end
