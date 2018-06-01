function [s] = collect4Dsmoothpow3(subjname)

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);
load([subjname,'powPrePreviousCong3']);
load('/home/jan/matlab/mri/templategrid6mm.mat');

addpath /home/jan/projects/ccc/3D/
for k = 1:numel(s)
  s{k}.stat1(isnan(s{k}.stat1)) = 0;
  s{k}.stat2(isnan(s{k}.stat2)) = 0;
  s{k}.stat3(isnan(s{k}.stat3)) = 0;
  s{k}.stat4(isnan(s{k}.stat4)) = 0;
end
for k = 1:numel(s{1}.freq)
  tmp = [];
  tmp.fwhm   = double(s{1}.fwhm(:,k));
  tmp.inside = double(s{1}.inside);
  tmp.dim    = double(s{1}.dim);
  tmp.pos    = double(grid.pos);
  krn = compute_kernel(tmp, 'truncate', 1e-6)';
  in  = s{1}.inside;
  for kk = 1:numel(s)
    nrm = (s{kk}.stat1+s{kk}.stat2+s{kk}.stat3+s{kk}.stat4)./4;
    s{kk}.stat1(in,k) = single(krn*double(s{kk}.stat1(in,k)./nrm(in,k)));
    s{kk}.stat2(in,k) = single(krn*double(s{kk}.stat2(in,k)./nrm(in,k)));
    s{kk}.stat3(in,k) = single(krn*double(s{kk}.stat3(in,k)./nrm(in,k)));
    s{kk}.stat4(in,k) = single(krn*double(s{kk}.stat4(in,k)./nrm(in,k)));
  end
  clear krn
end

save([subjname,'powPrePreviousCong3Smooth'], 's');
