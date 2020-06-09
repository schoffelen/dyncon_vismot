load list
for k=1:19
subject = vismot_subjinfo(list{k});

if isfield(subject, 'mriname') && ~isempty(subject.mriname)
    name = subject.mriname;
  else
    name = subject.name;
  end
  fprintf('processing subject %s\n',subject.name);

  cd(subject.rawpath);
  cd([name,'_RAW']);
  d   = dir([upper(name),'*']);
  mri = ft_read_mri(d(end).name);
  
  age{k} = mri.hdr(1).PatientAge;
  sex{k} = mri.hdr(1).PatientSex;
end
for k=1:19
a(k) = str2num(age{k}(2:3));
end