subjinfo


freqs   = [8:2:36 40:4:100];
datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/'; 
for k = 1:length(subjno)
for m = 1:length(freqno)
  
  subject       = SUBJ(subjno(k));
  frequency     = freqs(freqno(m));
  [statlc, statli, statrc, statri] = doSourceanalysisDICScohPre(subject, frequency);

  warning off;
  statlc = struct2single(statlc);
  statli = struct2single(statli);
  statrc = struct2single(statrc);
  statri = struct2single(statri);
  warning on;
  
  cd([datadir,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'cohPreCongClean.mat'],'statlc','statli','statrc','statri');
  clear statlc statli statrc statri;
end
end
