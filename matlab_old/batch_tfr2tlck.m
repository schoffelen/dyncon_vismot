subjinfo

sel = find(~ismember({SUBJ(:).name}, {'GAR12' 'BKA01' 'KBI24'}));

cd('/analyse/1/Project0002/tmpProject0030');
for k = 1:numel(sel)
  subject = SUBJ(sel(k))
  load([subject.name,'tfrperiHanning2fourier']);
  for m = 1:4
    [tlck, covariance{m}, dof{m}] = tfr2tlck(allfreq{m}, 'all');
  end
  freq = allfreq{1}.freq;
  save([subject.name,'tfr2tlck'], 'tlck', 'covariance', 'dof', 'freq');
  clear allfreq covariance dof tlck
end
