subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject         = SUBJ(subjno(k));
  frequency       = freqs(freqno(m));
  [stat,fwhm] = doSourceanalysisDICSprepst(subject,frequency);
  
  warning off;
  stat = struct2single(stat);
  warning on;
  
  %cd([subject.pathname,'source/']);
  cd('/analyse/1/Project0002/tmpProject0030/');
  save([subject.name,'stat',num2str(frequency,'%03d'),'prepst.mat'],'stat','fwhm');
  clear stat fwhm;
end
end
