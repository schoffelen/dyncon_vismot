subjinfo;
for k = [1:3 5:20]
  subject = SUBJ(k)
  if k>3,
    flag = 1;
  else
    flag = 0;
  end
  [output, input, binaxis] = stratifyRT(subject, flag);
  cd('/analyse/4/Project0030/stratifyRT');
  save([subject.name,'stratifyRT'],'output','input','binaxis');
  clear output input;
end
