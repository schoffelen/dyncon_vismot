subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject         = SUBJ(subjno(k));
  frequency       = freqs(freqno(m));
  [stat13,stat42,stat12,stat43,design] = doSourceanalysisDICSpst2(subject,frequency);
  
  warning off;
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  warning on;
  
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'pstPPI.mat'],'stat13','stat42','stat12','stat43','design');
  clear stat13 stat42 stat12 stat43 design;
end
end
