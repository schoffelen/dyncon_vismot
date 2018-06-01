subjinfo;
savepath = '/analyse/4/Project0030/freq/tfr';

band      = 'alpha';
tfwindow  = [128 128 128]./256;
frequency = [8 10 12];
smoothing = [0 0 0];
for k = 1:length(subjno)
  indx    = subjno(k);
  subject = SUBJ(indx);
  if indx<=3,
    flag = 0;
  else
    flag = 1;
  end
  freq    = doFreqanalysisTFR(subject,frequency,tfwindow,smoothing,flag);

  savename = [subject.name,'tfr',band];
  cd(savepath);
  save(savename, 'freq');
  clear freq;
end
