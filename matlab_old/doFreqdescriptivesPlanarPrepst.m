freqdir = '/analyse/4/Project0030/freq/';
cd(freqdir);
d     = dir;
d     = d(3:end);
names = {d(:).name}';
x     = strfind(names, 'peri');
sel   = find(~cellfun('isempty',x));
names = names(sel);
for k = 1:numel(names)
  names{k}
  load(names{k});
  for m = 1:4
    fd{k,m} = ft_freqdescriptives([], allfreq{m});
  end
  clear allfreq;
end
