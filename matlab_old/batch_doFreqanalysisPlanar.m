subjinfo;

for k = 1:length(subjno)
  subject = SUBJ(subjno(k));
  %[stat] = doFreqanalysisPlanar(subject,smoothing,foilim);
  [stat] = doFreqanalysisPlanar(subject,smoothing);

  cd(subject.pathname);
  cd('freq');
  save([subject.name,'planar',num2str(smoothing,'%03d')],'stat');
  %clear stat lrp allfreq;
  clear stat;
end
