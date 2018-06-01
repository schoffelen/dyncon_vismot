subjinfo;
for k = 1:numel(subjno)
  indx    = subjno(k);
  subject = SUBJ(indx);
  if indx>3,
    flag = 1;
  else
    flag = 0;
  end
  rt = data2rt(subject, flag);
  cd(subject.pathname);
  cd('rt');
  save([subject.name,'rt'], 'rt');
  clear rt;
end
