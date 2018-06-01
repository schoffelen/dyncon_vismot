subjinfo;

for k = 1:length(subjno)
  indx    = subjno(k);
  subject = SUBJ(indx);
  if indx>3,
    [freqpre] = doFreqanalysisMtmfftPre400(subject,smoothing,foilim,1);
    %[freqpre] = doFreqanalysisMtmfftPre625(subject,smoothing,foilim,1);
  else
    %data are not yet aligned
    [freqpre] = doFreqanalysisMtmfftPre400(subject,smoothing,foilim,0);
    %[freqpre] = doFreqanalysisMtmfftPre625(subject,smoothing,foilim,0);
  end

  for kk = 1:4
    warning off
    freqpre{kk} = struct2single(freqpre{kk});
    warning on
  end

  cd('/analyse/4/Project0030/');
  cd('freq');
  save([subject.name,'mtmfftPre',num2str(round(10*smoothing),'%03d'),'.mat'],'freqpre');
  %save([subject.name,'mtmfftPre625_',num2str(round(10*smoothing),'%03d'),'.mat'],'freqpre');
  clear freqpre;
end
