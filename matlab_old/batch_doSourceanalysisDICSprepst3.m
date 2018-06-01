subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject   = SUBJ(subjno(k));
  frequency = freqs(freqno(m));
  if subjno(k)>3
   [stat13,stat42,stat13b,stat42b] = doSourceanalysisDICSprepst3(subject,frequency,1);
  else
   [stat13,stat42,stat13b,stat42b] = doSourceanalysisDICSprepst3(subject,frequency);
  end
  warning off;
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  stat13b = struct2single(stat13b);
  stat42b = struct2single(stat42b);
  warning on;
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'prepstGLM.mat'],'stat13','stat42','stat13b','stat42b');
  %save([subject.name,'stat10hanning.mat'],'stat13','stat42');
  clear stat13 stat42 stat13b stat42b;
end
end
