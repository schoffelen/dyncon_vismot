subjinfo

freqs = [4:2:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  subject   = SUBJ(subjno(k));
  frequency = freqs(freqno(m));
  [stat13,stat42,stat13b,stat42b] = doSourceanalysisDICSglm4(subject,frequency);
  warning off;
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  stat13b = struct2single(stat13b);
  stat42b = struct2single(stat42b);
  warning on;
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'glm4.mat'],'stat13','stat42','stat13b','stat42b');
  clear stat13 stat42 stat13b stat42b;
end
end
