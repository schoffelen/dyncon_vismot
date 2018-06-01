subjinfo
for k = [1:3 5:20]
  subject = SUBJ(k);
  block   = event2block(subject);
  cd(subject.pathname);
  cd('block');
  save([subject.name,'block'], 'block');
  clear block
end
