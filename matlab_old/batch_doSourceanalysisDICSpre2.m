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

%freqs = [2:1:22 16:2:40 40:3:100];
freqs = [16:2:28];
for k = 1:length(subjno)
for m = 1:length(freqno)
  if freqno(m)<=0,
    smoothing = 0;
  elseif freqno(m)<=34,
    smoothing = 3.75;
    %smoothing = 5;
  else
    smoothing = 7.5;
  end
  
  subject       = SUBJ(subjno(k));
  frequency     = freqs(freqno(m));
  [statlc, statli, statrc, statri] = doSourceanalysisDICSpre2(subject, frequency, smoothing)


  warning off;
  statlc = struct2single(statlc);
  statli = struct2single(statli);
  statrc = struct2single(statrc);
  statri = struct2single(statri);
  warning on;
  
  cd([subject.pathname,'source/']);
  save([subject.name,'stat',num2str(frequency,'%03d'),'pre2roicong_',num2str(round(10*smoothing),'%03d'),'.mat'],'statlc','statli','statrc','statri');
  clear statlc statli statrc statri;
end
end
