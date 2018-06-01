subjinfo

for k = [1:3 5:20]
  subject = SUBJ(k);
  if k<=3
    [stat13,stat42,stat13b,stat42b] = doFreqanalysisPrePstStratified2(subject,[6 40], 0);
  else
    [stat13,stat42,stat13b,stat42b] = doFreqanalysisPrePstStratified2(subject,[6 40], 1);
  end
  savename = [subject.name,'glm004'];
  savepath = '/analyse/4/Project0030/freq/glm2';
  cd(savepath);
  save(savename,'stat13','stat42','stat13b','stat42b');
end


