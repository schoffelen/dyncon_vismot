function concatenatesourcePrePstGLM(subject,band,bandname)

cd([subject.pathname,'source/']);

basename = [subject.name, 'stat'];
d        = dir;
d        = d(3:end);
names    = {d.name};
sel      = strmatch(basename, names);
d        = d(sel);
names    = {d.name};
sel      = strfind(names, 'GLM');
sel      = find(cellfun(@isempty, sel)==0);
d        = d(sel);
names    = names(sel);
sel      = strfind(names, 'resplock');
sel      = find(cellfun(@isempty, sel)==0);
d        = d(sel);
if nargin>1,
  band = band(1):band(2);
  ok   = zeros(length(d),1);
  for k = 1:length(d)
    ok(k) = sum(band==str2num(d(k).name(end-23:end-21)));
  end
  d = d(ok==1);
end

%if strcmp(subject.name, 'JOE22') || strcmp(subject.name, 'MHH14') || strcmp(subject.name, 'T001')
%  d = d(1:2:end);
%else
%  d = d(1:3:end);
%end

for k = 1:length(d)
  fname = d(k).name;
  fprintf('loading %s\n',fname);
  load(fname);
  stat13.dim  = [size(stat13.pos,1) 1];
  stat42.dim  = [size(stat42.pos,1) 1];
  stat13b.dim = [size(stat13.pos,1) 1];
  stat42b.dim = [size(stat42.pos,1) 1];
  stat13.dimord = 'pos_freq';
  stat42.dimord = 'pos_freq';
  stat13b.dimord = 'pos_freq';
  stat42b.dimord = 'pos_freq';
  stat13.freq  = str2num(fname(end-23:end-21));
  stat42.freq  = str2num(fname(end-23:end-21));
  stat13b.freq = str2num(fname(end-23:end-21));
  stat42b.freq = str2num(fname(end-23:end-21));
  s13{k} = stat13;
  s42{k} = stat42;
  s13b{k} = stat13b;
  s42b{k} = stat42b;
  clear stat13 stat42 stat13b stat42b;
end
%stat    = selectdata(s{:},    'param', 'stat');
%stat13  = selectdata(s13{:},  'param', 'stat');
%stat42  = selectdata(s42{:},  'param', 'stat');
%stat13b = selectdata(s13b{:}, 'param', 'stat');
%stat42b = selectdata(s42b{:}, 'param', 'stat');
%stat13x = selectdata(s13x{:}, 'param', 'stat');
%stat42x = selectdata(s42x{:}, 'param', 'stat');

stat13  = selectdata(s13{:},  'param', {'stat2' 'stat2rt'});
stat42  = selectdata(s42{:},  'param', {'stat2' 'stat2rt'});
stat13b = selectdata(s13b{:}, 'param', {'stat2' 'stat2rt'});
stat42b = selectdata(s42b{:}, 'param', {'stat2' 'stat2rt'});

cd('4DprepstGLMresplock');
if nargin>1,
  fname = [subject.name,'stat4DprepstGLMresplock_',bandname];
else
  fname = [subject.name,'stat4Dprepst'];
end
save(fname, 'stat13', 'stat42', 'stat13b', 'stat42b');


