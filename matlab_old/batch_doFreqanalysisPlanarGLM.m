subjinfo;

for k = 1:length(subjno)
  subject = SUBJ(subjno(k));
  [stat]  = doFreqanalysisPlanarGLM(subject,smoothing,foilim);

  cd(subject.pathname);
  cd('freq');
  save([subject.name,'planarGLM',num2str(smoothing,'%03d')],'stat');
  clear stat;
end
