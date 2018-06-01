dir1 = '/analyse/4/Project0030/freq/glm';
dir2 = '/analyse/4/Project0030/freq/glm2';

cfg        = [];
cfg.layout = '4D248.lay';
lay        = prepare_layout(cfg);
label      = lay.label;

cd(dir1);
d = dir;
d = d(3:end);
for k = 1:19
  fname = d(k).name;
  load(fname);
  s13{k,1}  = stat13;
  s13b{k,1} = stat13b;
  s42{k,1}  = stat42;
  s42b{k,1} = stat42b;
end
cd(dir2);
d = dir;
d = d(3:end);
for k = 1:19
  fname = d(k).name;
  load(fname);
  s13{k,2}  = stat13;
  s13b{k,2} = stat13b;
  s42{k,2}  = stat42;
  s42b{k,2} = stat42b;
end

for k = 1:19
  x1(k,1) = mean(mean(s13{k}.stat));
  x2(k,1) = mean(mean(s42{k}.stat));
  x3(k,1) = mean(mean(s13b{k}.stat));
  x4(k,1) = mean(mean(s42b{k}.stat));
  x1(k,2) = mean(mean(s13{k}.stat2));
  x2(k,2) = mean(mean(s42{k}.stat2));
  x3(k,2) = mean(mean(s13b{k}.stat2));
  x4(k,2) = mean(mean(s42b{k}.stat2));
  x1(k,3) = mean(mean(s13{k}.stat3));
  x2(k,3) = mean(mean(s42{k}.stat3));
  x3(k,3) = mean(mean(s13b{k}.stat3));
  x4(k,3) = mean(mean(s42b{k}.stat3));
  x1b(k,1) = mean(mean(s13{k,2}.stat));
  x2b(k,1) = mean(mean(s42{k,2}.stat));
  x3b(k,1) = mean(mean(s13b{k,2}.stat));
  x4b(k,1) = mean(mean(s42b{k,2}.stat));
  x1b(k,2) = mean(mean(s13{k,2}.stat2));
  x2b(k,2) = mean(mean(s42{k,2}.stat2));
  x3b(k,2) = mean(mean(s13b{k,2}.stat2));
  x4b(k,2) = mean(mean(s42b{k,2}.stat2));
  x1b(k,3) = mean(mean(s13{k,2}.stat3));
  x2b(k,3) = mean(mean(s42{k,2}.stat3));
  x3b(k,3) = mean(mean(s13b{k,2}.stat3));
  x4b(k,3) = mean(mean(s42b{k,2}.stat3));
end
