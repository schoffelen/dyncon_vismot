subjinfo;
names = {SUBJ(:).name}';
sel   = find(~ismember(names, {'BKA01' 'GAR12' 'KBI24'}));

savepath = '/analyse/4/Project0030/stat';

for k = 1:numel(sel)
  subjno  = sel(k);
  subject = SUBJ(subjno);
  [stat13, stat42] = doFreqstatisticsPlanarTFRpre(subject);
  cd(savepath);
  save([subject.name,'planarTFRpre'], 'stat13', 'stat42');
  clear stat13 stat42;
end
