subjinfo


freqs   = [8:2:36 40:4:100];
datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/'; 
for k = 1:length(subjno)
for m = 1:length(freqno)
  
  subject       = SUBJ(subjno(k));
  frequency     = freqs(freqno(m));
  [statc, stati] = doSourceanalysisDICSpowPrePrevious2(subject, frequency);

  warning off;
  statc = struct2single(statc);
  stati = struct2single(stati);
  warning on;
  
  cd([datadir,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'powPrePreviousCong2.mat'],'statc','stati');
  clear statc stati;
end
end
