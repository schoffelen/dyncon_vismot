subjinfo;
for m = 1:length(subjno)
  subject = SUBJ(subjno(m));
  [data1,data2,data3,data4,data5] = doRealignment(subject);
  cd(subject.pathname);
  cd('data');
  save([subject.name,'data_aligned'],'data1','data2','data3','data4','data5');
  clear data1 data2 data3 data4 data5
end
