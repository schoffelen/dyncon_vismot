function [s13, s42, inside] = doSourcegrandavgCongSmooth(foi, band)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,'powCongSmooth')) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
    fprintf('loading %s\n', fname);
    load(fname);
    s{1}.avg.pow = s{1}.stat; s{1} = rmfield(s{1}, 'stat');
    s{2}.avg.pow = s{2}.stat; s{2} = rmfield(s{2}, 'stat');
    s{1} = ft_selectdata(s{1}, 'foilim', [foi(1) foi(end)]);
    s{2} = ft_selectdata(s{2}, 'foilim', [foi(1) foi(end)]);
    s{1}.pos = grid.pos
    s{2}.pos = grid.pos;
    cnt = cnt+1;
    s13{cnt} = s{1};
    s42{cnt} = s{2};
    if cnt == 1
      inside = zeros(s{1}.dim);
    end
    inside(s{1}.inside) = inside(s{1}.inside) + 1;
    names{end+1} = fname(1:5);
  end
end

s13 = ft_selectdata(s13{:}, 'param', 'pow');
s13.dim = s13.dim(end-2:end);
s42 = ft_selectdata(s42{:}, 'param', 'pow');
s42.dim = s42.dim(end-2:end);

save(['grandavgCongSmooth-',band], 's13', 's42', 'inside', 'names');
