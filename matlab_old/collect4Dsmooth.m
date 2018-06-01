function [s] = collect4Dsmooth(subjname)

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);
load([subjname,'acPreFastSlow']);
%load([subjname,'cohFastSlow']);

for k = 1:4
  s{k}.coh1(isnan(s{k}.coh1)) = 0;
  s{k}.coh2(isnan(s{k}.coh2)) = 0;
  s{k}.coh3(isnan(s{k}.coh3)) = 0;
  s{k}.coh4(isnan(s{k}.coh4)) = 0;
end
for k = 1:numel(s{1}.freq)
  tmp = [];
  tmp.fwhm   = double(s{1}.fwhm(:,k));
  tmp.inside = double(s{1}.inside);
  tmp.dim    = double(s{1}.dim);
  tmp.pos    = double(s{1}.pos);
  krn = compute_kernel(tmp, 'truncate', 1e-6)';
  in  = s{1}.inside;
  for kk = 1:4
    s{kk}.coh1(in,k) = single(krn*double(s{kk}.coh1(in,k)));
    s{kk}.coh2(in,k) = single(krn*double(s{kk}.coh2(in,k)));
    s{kk}.coh3(in,k) = single(krn*double(s{kk}.coh3(in,k)));
    s{kk}.coh4(in,k) = single(krn*double(s{kk}.coh4(in,k)));
  end
  clear krn
end

%save([subjname,'dcohFastSlowSmooth'], 's');
save([subjname,'cohFastSlowSmooth'], 's');
