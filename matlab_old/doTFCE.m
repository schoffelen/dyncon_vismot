function doTFCE(pathname, suffix, saveflag)

%only works on field stat2, and if input file contains variables stat13 and stat42
cd(pathname);
d     = dir;
names = {d(3:end).name}';
sel   = ~cellfun(@isempty, strfind(names, suffix));
names = names(sel);
for k = 1:numel(names)
  names{k}
  load(names{k});
  dim = pos2dim3d(stat13.pos);
  h0  = median(stat13.stat2(stat13.inside,:));
  tmp = stat13.stat2;
  tmp(stat13.inside,:) = h0(ones(numel(stat13.inside),1),:);
  tmp = reshape(tmp, [dim numel(stat13.freq)]);
  stat13.stat2 = tfce(reshape(stat13.stat2, [dim numel(stat13.freq)]),[],[],[],tmp);
  h0  = median(stat42.stat2(stat42.inside,:));
  tmp = stat42.stat2;
  tmp(stat42.inside,:) = h0(ones(numel(stat42.inside),1),:);
  tmp = reshape(tmp, [dim numel(stat42.freq)]);
  stat42.stat2 = tfce(reshape(stat42.stat2, [dim numel(stat42.freq)]),[],[],[],tmp);
  if saveflag,
    %save([names{k},'tfce'], 'stat13', 'stat42');
    save([names{k},'tfce2'], 'stat13', 'stat42');
  end
  clear stat13 stat42
end
