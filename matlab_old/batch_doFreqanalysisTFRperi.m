subjinfo

%frequency = [10:1:25];
%twindow   = 4./frequency;
%smoothing = zeros(1,numel(frequency));
frequency = [8:2:24];
twindow   = ones(1,numel(frequency))*0.5;
smoothing = zeros(1,numel(frequency));
toi       = [-192:4:192]./256;

for k = 1:numel(subjno)
  subject = SUBJ(subjno(k));
  [allfreq] = doFreqanalysisTFRperi(subject,frequency,smoothing,twindow,toi);
  %cd(subject.pathname);
  %cd('freq');
  cd /analyse/1/Project0002/tmpProject0030
  %save([subject.name,'tfrperiHanning'],'allfreq');
  save([subject.name,'tfrperiHanning2fourier'],'allfreq');
  clear allfreq;
end
out = 1;
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
