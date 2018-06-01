function [stat] = processLCMV(subject)

%this function loads the output of doSourceanalysisLCMV and processes the
%'pos_time' matrices: non-homogeneous smoothing + pooling across contrasts,
%where contrast for motor response left is flipped along the midline

cd([subject.pathname,'source_lcmv']);
load([subject.name,'lcmv'], 'stat');

if isfield(stat, 'statx')
  %nothing to do
else
  %convert to double precision
  warning off;
  stat = struct2double(stat);
  warning on;

  %compute smoothing kernel
  addpath /home/jan/projects/ccc/3D
  fprintf('computing smoothing kernel\n');
  insideold   = stat.inside;
  insidenew   = find(isfinite(stat.fwhm));
  stat.inside = insidenew; 
  krn         = compute_kernel(stat, 'truncate', 2.5e-5, 'feedback', 1);

  %take absolute value after baseline subtraction
  fprintf('preprocessing data\n');
  sel         = nearest(stat.time, 0);
  stat.stat13 = abs(blc(stat.stat13, [1 sel]));
  stat.stat42 = abs(blc(stat.stat42, [1 sel]));

  %smooth data
  fprintf('smoothing data\n');
  dim         = pos2dim3d(stat.pos);
  stat.stat13 = smooth_vol(stat.stat13, krn, dim, insidenew);
  stat.stat42 = smooth_vol(stat.stat42, krn, dim, insidenew);

  %pool data
  fprintf('pooling data\n');
  stat.statx = zeros(size(stat.stat13));
  for k = 1:length(stat.time)
    dum1 = reshape(stat.stat13(:,k), dim);
    dum2 = reshape(stat.stat42(:,k), dim);
    if all(dum1(:)==0) || any(~isfinite(dum1(:))),
      stat.statx(:,k) = reshape(dum2, prod(dim), 1);
    elseif all(dum2(:)==0) || any(~isfinite(dum2(:))),
      stat.statx(:,k) = reshape(flipdim(dum1,1), prod(dim), 1);
    else
      stat.statx(:,k) = (reshape(dum2, prod(dim), 1) + ...
                         reshape(flipdim(dum1,1), prod(dim), 1))./sqrt(2);
    end
  end

  warning off;
  stat = struct2single(stat);
  warning on;
  save([subject.name,'lcmv'], 'stat');
end
