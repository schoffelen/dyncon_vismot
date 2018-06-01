subjinfo

%savedir = '/analyse/4/Project0030/freqstats';
savedir = '/analyse/1/Project0002/tmpProject0030/';
names   = {SUBJ(:).name}';
sel     = find(~ismember(names, {'GAR12' 'BKA01' 'KBI24'}));
for k = 1:numel(sel)
%for k = 15:numel(sel)
  subject   = SUBJ(sel(k));
  bas       = [-0.3 -0.25];
  [avg,cnd] = doFreqstatsPlanarTFRperiRelchange(subject,bas);
  cd(savedir);
  save([subject.name,'relchange'], 'avg', 'cnd');
  clear avg cnd
end
