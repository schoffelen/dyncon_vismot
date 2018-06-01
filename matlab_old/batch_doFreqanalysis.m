subjinfo;

for k = 1:length(subjno)
  subject = SUBJ(subjno(k));
  [freqpre, freqpst] = doFreqanalysis(subject,smoothing,foilim);

  cd(subject.pathname);
  cd('freq');
  save([subject.name,'mtmfft',num2str(smoothing,'%03d')],'freqpre','freqpst'); %tapsmo=4
  clear freqpre freqpst;
end
