function [s] = collect4DsmoothpowPst(subjname)

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);
load([subjname,'powPstCong']);
%load([subjname,'powCong']);

for k = 1:2
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
  for kk = 1:2
    s{kk}.stat(in,k) = single(krn*double(s{kk}.stat(in,k)));
  end
  clear krn
end

save([subjname,'powCongSmooth'], 's');
