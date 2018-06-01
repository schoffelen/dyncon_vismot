subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject       = SUBJ(subjno(k));
  frequency     = freqs(freqno(m));
  [stat,design] = doSourceanalysisDICSpst3(subject,frequency);
  
  warning off;
  stat = struct2single(stat);
  warning on;
  
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'pstPPI2.mat'],'stat','design');
  clear stat design;
end
end
