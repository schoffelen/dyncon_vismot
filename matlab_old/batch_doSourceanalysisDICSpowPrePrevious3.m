subjinfo


%freqs   = [8:2:36 40:4:100];
freqs   = [4:2:14];
datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/'; 
for k = 1:length(subjno)
for m = 1:length(freqno)
  
  subject       = SUBJ(subjno(k));
  frequency     = freqs(freqno(m));
  [statl, statr] = doSourceanalysisDICSpowPrePrevious3(subject, frequency);

  warning off;
  statl = struct2single(statl);
  statr = struct2single(statr);
  warning on;
  
  cd([datadir,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'powPrePreviousCong3.mat'],'statl','statr');
  clear statl statr;
end
end
