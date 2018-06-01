function [stat13, stat42] = concatenatesource(subject)

cd([subject.pathname,'source/']);

basename = [subject.name, 'stat'];
d        = dir;
d        = d(3:end);
names    = {d.name};
sel      = strmatch(basename, names);
d        = d(sel);
if strcmp(subject.name, 'JOE22') || strcmp(subject.name, 'MHH14') || strcmp(subject.name, 'T001')
  d = d(1:2:end);
else
  d = d(1:3:end);
end

for k = 1:length(d)
  fname = d(k).name;
  fprintf('loading %s\n',fname);
  load(fname);
  stat13.dim = [size(stat13.pos,1) 1];
  stat42.dim = [size(stat42.pos,1) 1];
  stat13.dimord = 'pos_freq';
  stat42.dimord = 'pos_freq';
  %stat13.freq = str2num(fname(end-9:end-7));
  %stat42.freq = str2num(fname(end-9:end-7));
  %stat13.freq = str2num(fname(end-13:end-11));
  %stat42.freq = str2num(fname(end-13:end-11));
  stat13.freq = str2num(fname(end-6:end-4));
  stat42.freq = str2num(fname(end-6:end-4));
  s13{k} = stat13;
  s42{k} = stat42;
  clear stat13 stat42;
end
stat13 = selectdata(s13{:}, 'param', 'stat');
stat42 = selectdata(s42{:}, 'param', 'stat');

cd('4D');
fname = [subject.name,'stat4DRough'];
save(fname, 'stat13', 'stat42');


