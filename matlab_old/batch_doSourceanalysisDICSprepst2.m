subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject   = SUBJ(subjno(k));
  frequency = freqs(freqno(m));
  if subjno(k)>3
   [stat,stat13,stat42,stat13b,stat42b,stat13x,stat42x] = doSourceanalysisDICSprepst2(subject,frequency,1);
  else
   [stat,stat13,stat42,stat13b,stat42b,stat13x,stat42x] = doSourceanalysisDICSprepst2(subject,frequency);
  end
  warning off;
  stat   = struct2single(stat);
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  stat13b = struct2single(stat13b);
  stat42b = struct2single(stat42b);
  stat13x = struct2single(stat13x);
  stat42x = struct2single(stat42x);
  warning on;
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'prepst.mat'],'stat','stat13','stat42','stat13b','stat42b','stat13x','stat42x');
  %save([subject.name,'stat10hanning.mat'],'stat13','stat42');
  clear stat13 stat42 stat13b stat42b stat13x stat42x stat;
end
end
