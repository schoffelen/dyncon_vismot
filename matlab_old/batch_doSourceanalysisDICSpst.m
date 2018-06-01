subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject         = SUBJ(subjno(k));
  frequency       = freqs(freqno(m));
  [stat13,stat42] = doSourceanalysisDICSpst(subject,frequency);
  
  warning off;
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  warning on;
  
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'pst.mat'],'stat13','stat42');
  clear stat13 stat42 stat13b stat42b;
end
end
