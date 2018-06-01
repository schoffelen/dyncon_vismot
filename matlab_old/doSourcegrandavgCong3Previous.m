function [sl, sr, inside] = doSourcegrandavgCong3Previous(foi, band)

%cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/source/');
d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,'powPrePreviousCong3Smooth')) && isempty(strfind(d(k).name,'stat')) && isempty(strfind(d(k).name,'grandavg')),
  %if ~isempty(strfind(d(k).name,'powPrePreviousCongHann')) && isempty(strfind(d(k).name,'stat')) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
    fprintf('loading %s\n', fname);
    load(fname);
    s{1}.avg.pow(:,:,1) = s{1}.stat1; s{1} = rmfield(s{1}, 'stat1');
    s{1}.avg.pow(:,:,2) = s{1}.stat2; s{1} = rmfield(s{1}, 'stat2');
    s{1}.avg.pow(:,:,3) = s{1}.stat3; s{1} = rmfield(s{1}, 'stat3');
    s{1}.avg.pow(:,:,4) = s{1}.stat4; s{1} = rmfield(s{1}, 'stat4');
    s{2}.avg.pow(:,:,1) = s{2}.stat1; s{2} = rmfield(s{2}, 'stat1');
    s{2}.avg.pow(:,:,2) = s{2}.stat2; s{2} = rmfield(s{2}, 'stat2');
    s{2}.avg.pow(:,:,3) = s{2}.stat3; s{2} = rmfield(s{2}, 'stat3');
    s{2}.avg.pow(:,:,4) = s{2}.stat4; s{2} = rmfield(s{2}, 'stat4');
    s{1}.pos = grid.pos;
    s{2}.pos = grid.pos;
    s{1}.dimord = 'pos_freq_time';
    s{2}.dimord = 'pos_freq_time';
    s{1}.time = [1 2 3 4];
    s{2}.time = [1 2 3 4];
    s{1} = ft_selectdata(s{1}, 'foilim', [foi(1) foi(end)]);
    s{2} = ft_selectdata(s{2}, 'foilim', [foi(1) foi(end)]);
    cnt = cnt+1;
    sl{cnt} = s{1};
    sr{cnt} = s{2};
    if cnt == 1
      inside = zeros(s{1}.dim);
    end
    inside(s{1}.inside) = inside(s{1}.inside) + 1;
    names{end+1} = fname(1:5);
  end
end

sl = ft_selectdata(sl{:}, 'param', 'pow');
sl.dim = sl.dim(end-2:end);
sr = ft_selectdata(sr{:}, 'param', 'pow');
sr.dim = sr.dim(end-2:end);

save(['grandavgPrePreviousCong3Smooth-',band], 'sl', 'sr', 'inside', 'names');
