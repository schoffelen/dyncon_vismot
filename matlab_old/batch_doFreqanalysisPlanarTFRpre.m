subjinfo;

frequency = [10:1:25];
twindow   = 4./frequency;
smoothing = zeros(1,numel(frequency));

for k = 1:numel(subjno)
  subject = SUBJ(subjno(k));
  [allfreq] = doFreqanalysisPlanarTFRpre(subject,frequency,smoothing,twindow);

  cd(subject.pathname);
  cd('freq');
  save([subject.name,'tfrpreHanning'],'allfreq');
  clear allfreq;
end

%frequency = 20;
%smoothing = 4;
%twindow   = 0.5;
%
%for k = 1:length(subjno)
%  subject = SUBJ(subjno(k));
%  [allfreq] = doFreqanalysisPlanarTFR(subject,frequency,smoothing,twindow);
%
%  cd(subject.pathname);
%  cd('freq');
%  save([subject.name,'tfr',num2str(frequency,'%03d')],'allfreq');
%  clear allfreq;
%end

%frequency = 56;
%smoothing = 12;
%twindow   = 0.250;
%
%for k = 1:length(subjno)
%  subject = SUBJ(subjno(k));
%  [allfreq] = doFreqanalysisPlanarTFR(subject,frequency,smoothing,twindow);
%
%  cd(subject.pathname);
%  cd('freq');
%  save([subject.name,'tfr',num2str(frequency,'%03d')],'allfreq');
%  clear allfreq;
%end
