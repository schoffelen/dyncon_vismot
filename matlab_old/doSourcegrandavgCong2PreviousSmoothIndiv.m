function [sc, si, inside] = doSourcegrandavgCong2PreviousSmoothIndiv(foi, band)

%cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,'powPrePreviousCong4Smooth')) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
    fprintf('loading %s\n', fname);
    load(fname);
    cnt = cnt+1;
    s{1}.avg.pow = s{1}.stat; s{1} = rmfield(s{1}, 'stat');
    s{2}.avg.pow = s{2}.stat; s{2} = rmfield(s{2}, 'stat');
    s{1} = ft_selectdata(s{1}, 'foilim', [foi{1}(cnt) foi{1}(cnt)]);
    s{2} = ft_selectdata(s{2}, 'foilim', [foi{2}(cnt) foi{2}(cnt)]);
    s{1}.pos = grid.pos
    s{2}.pos = grid.pos;
    s{1}.freq = mean(foi{1});
    s{2}.freq = mean(foi{2});
    sc{cnt} = s{1};
    si{cnt} = s{2};
    if cnt == 1
      inside = zeros(s{1}.dim);
    end
    inside(s{1}.inside) = inside(s{1}.inside) + 1;
    names{end+1} = fname(1:5);
  end
end

sc = ft_selectdata(sc{:}, 'param', 'pow');
sc.dim = sc.dim(end-2:end);
si = ft_selectdata(si{:}, 'param', 'pow');
si.dim = si.dim(end-2:end);

save(['grandavgPrePreviousCong4Smooth-',band], 'sc', 'si', 'inside', 'names');
