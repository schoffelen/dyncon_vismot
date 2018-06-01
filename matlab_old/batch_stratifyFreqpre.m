subjinfo;

for m = 1:length(subjno)
  subject = SUBJ(subjno(m));
  fname   = [subject.pathname,'freq/',subject.name,'mtmfft',num2str(smoothing,'%03d')];
  load(fname, 'freqpre');
  [oktrials, label, freqs] = stratifyFreqpre(freqpre);
  fname    = [subject.pathname,'freq/',subject.name,'stratify',num2str(smoothing,'%03d')];
  save(fname,'oktrials','label','freqs');
end
