pathname = '/analyse/4/Project0030/source/';
foi      = 8:2:40;
for k = 1:numel(foi)
  foi(k)
  suffix = ['stat',num2str(foi(k),'%03d'),'pre2pow_038'];
  doplot = 0;
  nrand  = 1000;
  fname  = {'stat2'};
  collectVoxelstats2(pathname,suffix,doplot,nrand,fname);
end

cd(pathname);
d     = dir;
names = {d(3:end).name}';
sel   = strfind(names, 'grandavg');
sel   = ~cellfun('isempty', sel);
names = names(sel);
sel   = strfind(names, 'pre2pow');
sel   = ~cellfun('isempty', sel);
names = names(sel);
for k = 1:numel(names)
  load(names{k}, 'allstat');
  for m = 1:8
    tmp{m}{k} = allstat(m);
    tmp{m}{k}.avg.pow = tmp{m}{k}.stat;
    tmp{m}{k}.avg.nai = tmp{m}{k}.prob;
    tmp{m}{k} = rmfield(tmp{m}{k}, 'stat');
    tmp{m}{k} = rmfield(tmp{m}{k}, 'prob');
  end
end

for k = 1:8
  stat{k} = selectdata(tmp{k}{:}, 'param', 'pow');
  stat{k}.dimord = 'pos_freq';
  stat{k}.pow    = stat{k}.pow';
  stat{k}.freq   = foi;
end
for k = 1:8
  stat2{k} = selectdata(tmp{k}{:}, 'param', 'nai');
  stat2{k}.dimord = 'pos_freq';
  stat2{k}.nai    = stat2{k}.nai';
  stat2{k}.freq   = foi;
end

cd(pathname);
d     = dir;
names = {d(3:end).name}';
sel   = strfind(names, 'grandavg');
sel   = ~cellfun('isempty', sel);
names = names(sel);
sel   = strfind(names, 'pre2pow');
sel   = ~cellfun('isempty', sel);
names = names(sel);
for k = 1:numel(names)
  load(names{k}, 'stat');
  str = {'statlcstat2','statlistat2','statrcstat2','statristat2'};
  for m = 1:4
    tmp{m}{k} = stat{m};
    tmp{m}{k}.pow = tmp{m}{k}.(str{m});
    tmp{m}{k} = rmfield(tmp{m}{k}, str{m});
    tmp{m}{k}.freq = foi(k);
    tmp{m}{k}.powdimord = [stat{m}.([str{m},'dimord']),'_freq'];
    tmp{m}{k} = rmfield(tmp{m}{k}, [str{m},'dimord']);
  end
end
stat      = selectdata(tmp{4}{:}, 'param', 'pow');
stat.freq = foi;
stat.dim  = pos2dim3d(stat.pos);


