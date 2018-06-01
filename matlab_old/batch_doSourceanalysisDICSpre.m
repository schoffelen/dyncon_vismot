subjinfo

%smoothing = 0;
%if smoothing==0,
%  freqs = [2:1:25];
%end

%smoothing = 3.75;
%if smoothing==3.75,
%  freqs = [16:2:40];
%end

%smoothing = 7.5;
%if smoothing==7.5,
%  freqs = [40:3:100];
%end

freqs = [2:1:22 16:2:40 40:3:100];
for k = 1:length(subjno)
for m = 1:length(freqno)
  if freqno(m)<=21,
    smoothing = 0;
  elseif freqno(m)<=34,
    smoothing = 3.75;
  else
    smoothing = 7.5;
  end
  
  subject         = SUBJ(subjno(k));
  frequency       = freqs(freqno(m));
  [stat13,stat42] = doSourceanalysisDICSpre(subject,frequency,smoothing);
  
  warning off;
  stat13 = struct2single(stat13);
  stat42 = struct2single(stat42);
  warning on;
  
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'pre',num2str(round(10*smoothing),'%03d'),'.mat'],'stat13','stat42');
  clear stat13 stat42;
end
end
