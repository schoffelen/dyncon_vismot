subjinfo;

for k = 1:length(subjno)
  subject = SUBJ(subjno(k));
  if k>3,
    [freqpre, freqpst] = doFreqanalysisBalancePrePost2(subject,smoothing,foilim,1);
  else
    %data are not yet aligned
    [freqpre, freqpst] = doFreqanalysisBalancePrePost2(subject,smoothing,foilim,0);
  end

  for k = 1:5
    warning off
    freqpre{k} = struct2single(freqpre{k});
    freqpst{k} = struct2single(freqpst{k});
    warning on
  end

  cd(subject.pathname);
  cd('freq');
  if k>3,
    save([subject.name,'mtmfft_aligned',num2str(smoothing,'%03d'),'bnew.mat'],'freqpre','freqpst'); %tapsmo=4
  else
    save([subject.name,'mtmfft',num2str(smoothing,'%03d'),'bnew.mat'],'freqpre','freqpst'); %tapsmo=4
  end
  clear freqpre freqpst;
end
