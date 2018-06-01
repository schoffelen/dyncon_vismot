function [allroi,allfoi,names] = plotROI(suffix);

cd('~/matlab/mri/');
load('templategrid6mm');
cd('/analyse/4/Project0030/roi');
d     = dir;
names = {d(3:end).name}';
sel   = find(~cellfun('isempty', strfind(names, suffix)));
names = names(sel);

allroi = zeros([grid.dim numel(sel)]);
tmp    = zeros(grid.dim);
for k = 1:numel(sel)
  tmp(:) = 0;
  load(names{k});
  if numel(foi)==1,
    foi = [foi foi];
  end
  tmp(roi1) = 1;
  tmp(roi2) = 1;
  allroi(:,:,:,k) = tmp;
  %allfoi(k,:)       = foi;
  allfoi{k}       = foi;
end
