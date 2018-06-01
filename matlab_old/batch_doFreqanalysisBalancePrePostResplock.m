subjinfo;

for k = 1:length(subjno)
  subject = SUBJ(subjno(k));
  if subjno(k)>3,
    [freqpre, freqpst] = doFreqanalysisBalancePrePostResplock(subject,smoothing,foilim,1);
  else
    %data are not yet aligned
    [freqpre, freqpst] = doFreqanalysisBalancePrePostResplock(subject,smoothing,foilim,0);
  end

  for kk = 1:length(freqpre)
    warning off
    freqpre{kk} = struct2single(freqpre{kk});
    freqpst{kk} = struct2single(freqpst{kk});
    warning on
  end

  cd(subject.pathname);
  cd('freq');
  if subjno(k)>3,
    save([subject.name,'mtmfft_aligned',num2str(smoothing,'%03d'),'resplock.mat'],'freqpre','freqpst'); %tapsmo=4
  else
    save([subject.name,'mtmfft',num2str(smoothing,'%03d'),'resplock.mat'],'freqpre','freqpst'); %tapsmo=4
  end
  clear freqpre freqpst;
end
