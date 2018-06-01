function [s] = collect4Dsmoothpow(subjname)

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);
load([subjname,'powPrePreviousCong4']);

addpath /home/jan/projects/ccc/3D/
for k = 1:numel(s)
  s{k}.stat(isnan(s{k}.stat)) = 0;
end
for k = 1:numel(s{1}.freq)
  tmp = [];
  tmp.fwhm   = double(s{1}.fwhm(:,k));
  tmp.inside = double(s{1}.inside);
  tmp.dim    = double(s{1}.dim);
  tmp.pos    = double(s{1}.pos);
  krn = compute_kernel(tmp, 'truncate', 1e-6)';
  in  = s{1}.inside;
  for kk = 1:numel(s)
    s{kk}.stat(in,k) = single(krn*double(s{kk}.stat(in,k)));
  end
  clear krn
end

%save([subjname,'dcohFastSlowSmooth'], 's');
save([subjname,'powPrePreviousCong4Smooth'], 's');
