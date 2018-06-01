subjinfo

%savedir = '/analyse/4/Project0030/freqstats';
savedir = '/analyse/1/Project0002/tmpProject0030/';
names   = {SUBJ(:).name}';
sel     = find(~ismember(names, {'GAR12' 'BKA01' 'KBI24'}));
for k = 1:numel(sel)
%for k = 15:numel(sel)
  bas = [-0.3 -0.25];
  subject   = SUBJ(sel(k));
  [stat13,stat42,stat,stat13strat,stat42strat,statstrat] = doFreqstatsPlanarTFRperiRelCong(subject, bas);
  cd(savedir);
  save([subject.name,'relcong'], 'stat13', 'stat42', 'stat', 'stat13strat', 'stat42strat', 'statstrat');
  clear stat13 stat42 stat stat13strat stat42strat statstrat
end
