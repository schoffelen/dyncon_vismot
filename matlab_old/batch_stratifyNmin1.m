subjinfo

savepath = '/analyse/4/Project0030/dataPreStratify';

for k = [1:3 5:20]
  subject = SUBJ(k);
  if k<4,
    flag = 0;
  else
    flag = 1;
  end
  data = stratifyNmin1(subject,flag);
  cd(savepath);
  save([subject.name,'stratifyNmin1'], 'data');
  clear data
end
