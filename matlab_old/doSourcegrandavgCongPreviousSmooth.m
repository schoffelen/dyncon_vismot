function [slc, sli, src, sri, inside] = doSourcegrandavgCongPreviousSmooth(foi, band)

cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,'powPrePreviousCongSmooth')) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
    fprintf('loading %s\n', fname);
    load(fname);
    s{1}.avg.pow = s{1}.stat; s{1} = rmfield(s{1}, 'stat');
    s{2}.avg.pow = s{2}.stat; s{2} = rmfield(s{2}, 'stat');
    s{3}.avg.pow = s{3}.stat; s{3} = rmfield(s{3}, 'stat');
    s{4}.avg.pow = s{4}.stat; s{4} = rmfield(s{4}, 'stat');
    s{1} = ft_selectdata(s{1}, 'foilim', [foi(1) foi(end)]);
    s{2} = ft_selectdata(s{2}, 'foilim', [foi(1) foi(end)]);
    s{3} = ft_selectdata(s{3}, 'foilim', [foi(1) foi(end)]);
    s{4} = ft_selectdata(s{4}, 'foilim', [foi(1) foi(end)]);
    s{1}.pos = grid.pos
    s{2}.pos = grid.pos;
    s{3}.pos = grid.pos;
    s{4}.pos = grid.pos;
    cnt = cnt+1;
    slc{cnt} = s{1};
    sli{cnt} = s{2};
    src{cnt} = s{3};
    sri{cnt} = s{4};
    if cnt == 1
      inside = zeros(s{1}.dim);
    end
    inside(s{1}.inside) = inside(s{1}.inside) + 1;
    names{end+1} = fname(1:5);
  end
end

slc = ft_selectdata(slc{:}, 'param', 'pow');
slc.dim = slc.dim(end-2:end);
sli = ft_selectdata(sli{:}, 'param', 'pow');
sli.dim = sli.dim(end-2:end);
src = ft_selectdata(src{:}, 'param', 'pow');
src.dim = src.dim(end-2:end);
sri = ft_selectdata(sri{:}, 'param', 'pow');
sri.dim = sri.dim(end-2:end);

save(['grandavgPrePreviousCongSmooth-',band], 'slc', 'sli', 'src', 'sri', 'inside', 'names');
