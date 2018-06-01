function concatenatesourcePre(subject,band,bandname,smoothing)

% concatenatesourcePre(subject,band,bandname,smoothing)

cd([subject.pathname,'source/']);

basename = [subject.name, 'stat'];
d        = dir;
d        = d(3:end);
names    = {d.name};
sel      = strmatch(basename, names);
d        = d(sel);
names    = {d.name};
sel      = strfind(names, 'pre');
sel      = find(cellfun(@isempty, sel)==0);
d        = d(sel);
names    = names(sel);
sel      = strfind(names, smoothing);
sel      = find(cellfun(@isempty, sel)==0);
d        = d(sel);
names    = names(sel);
sel      = strfind(names, 'pst');
sel      = find(cellfun(@isempty, sel)==1);
d        = d(sel);
if nargin>1,
  band = band(1):band(2);
  ok   = zeros(length(d),1);
  for k = 1:length(d)
    ok(k) = sum(band==str2num(d(k).name(end-12:end-10)));
  end
  d = d(ok==1);
end

for k = 1:length(d)
  fname = d(k).name;
  fprintf('loading %s\n',fname);
  load(fname);
  stat13.dim  = [size(stat13.pos,1) 1];
  stat42.dim  = [size(stat42.pos,1) 1];
  stat13.dimord = 'pos_freq';
  stat42.dimord = 'pos_freq';
  stat13.freq  = str2num(fname(end-12:end-10));
  stat42.freq  = str2num(fname(end-12:end-10));
  s13{k} = stat13;
  s42{k} = stat42;
  clear stat13 stat42;
end
%stat13  = selectdata(s13{:},  'param', 'stat');
%stat42  = selectdata(s42{:},  'param', 'stat');

stat13  = selectdata(s13{:},  'param', {'stat2' 'stat2y'});
stat42  = selectdata(s42{:},  'param', {'stat2' 'stat2y'});

cd('4Dpre');
fname = [subject.name,'stat4Dpre',bandname];
save(fname, 'stat13', 'stat42');


