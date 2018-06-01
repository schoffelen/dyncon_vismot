subjinfo

for k = 1:length(subjno)
  subject   = SUBJ(subjno(k));
  if subjno(k)<=3,
    flag = 0;
  else
    flag = 1;
  end
  stat = doSourceanalysisLCMV(subject,flag);
  warning off;
  stat = struct2single(stat);
  warning on;
  cd([subject.pathname,'source_lcmv/']);
  save([subject.name,'lcmv.mat'],'stat');
  clear stat;
end
