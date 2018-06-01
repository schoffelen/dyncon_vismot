subjinfo;

for k = 1:length(subjno)
  indx    = subjno(k);
  subject = SUBJ(indx);
  if k>indx,
    [freqpre, freqpst] = doFreqanalysisMtmfft(subject,smoothing,foilim,1);
  else
    %data are not yet aligned
    [freqpre, freqpst] = doFreqanalysisMtmfft(subject,smoothing,foilim,0);
  end

  for kk = 1:4
    warning off
    freqpre{kk} = struct2single(freqpre{kk});
    freqpst{kk} = struct2single(freqpst{kk});
    warning on
  end

  cd(subject.pathname);
  cd('freq');
  save([subject.name,'mtmfft',num2str(smoothing,'%03d'),'.mat'],'freqpre','freqpst');
  clear freqpre freqpst;
end
