function concatenatesourcePrePst(subject,band,bandname)

cd([subject.pathname,'source/']);

basename = [subject.name, 'stat'];
d        = dir;
d        = d(3:end);
names    = {d.name};
sel      = strmatch(basename, names);
d        = d(sel);
names    = {d.name};
sel      = strfind(names, 'prepst');
sel      = find(cellfun(@isempty, sel)==0);
d        = d(sel);
names    = names(sel);
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
  stat.dim    = [size(stat.pos,1) 1];
  stat.dimord = 'pos_freq';
  stat.freq   = str2num(fname(end-12:end-10));
  stat.fwhm   = fwhm(:);
  s{k}        = stat;
  clear stat;
end

stat = selectdata(s{:},  'param', {'stat2' 'stat2y' 'fwhm'});

cd('4Dprepst');
fname = [subject.name,'stat4Dprepst',bandname];
save(fname, 'stat');


