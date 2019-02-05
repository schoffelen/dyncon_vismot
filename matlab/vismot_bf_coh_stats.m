for k = 1:numel(fnames)
  load(fnames{k});
  Z13(:,:,k) = zx13;
  Z42(:,:,k) = zx42;
  k
end

n = numel(fnames);
design = [ones(1,n) ones(1,n)*2;1:n 1:n];


cfg2.dim    = sourcemodel.dim;
cfg2.inside = find(sourcemodel.inside);
cfg2.threshold = 3.5;

dum = zeros(prod(sourcemodel.dim));
inside = find(sourcemodel.inside);
for k = 1:size(Z42,3)
  dum(inside,inside)=Z42(:,:,k);
  dum = reshape(dum,[sourcemodel.dim sourcemodel.dim]);
  dum = flip(dum,1);dum = flip(dum,4);
  dum = reshape(dum, [prod(sourcemodel.dim) prod(sourcemodel.dim)]);
  Z42(:,:,k) = dum(inside,inside);
end

%dum = Z13-nanmean(Z13,1);
%dum = dum-nanmean(Z13,2);
%dat = [reshape(dum,[],19) zeros(size(Z13,1)^2,19)];
%dat = [reshape(Z13,[],19) zeros(size(Z13,1)^2,19)];
dat = [reshape((Z13+Z42)./2,[],19) zeros(size(Z13,1)^2,19)];



indx = reshape(1:(size(Z13,1)^2),[],size(Z13,1));
indx = tril(indx,-1);
indx = indx(indx>0);
dat  = dat(indx,:);

cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
stat = ft_statfun_depsamplesT(cfg,dat,design);
stat.stat = squareform(stat.stat);
C         = vismot_cluster6D(cfg2,stat.stat);

curr_dir = pwd;
cd ~/matlab/fieldtrip/private
cfg3 = [];
cfg3.numrandomization = 500;
cfg3.ivar = 1;
cfg3.uvar = 2;
cfg3.resampling = 'permutation';
permutations = resampledesign(cfg3, design);
cd(curr_dir);

for k = 1:size(permutations,1)
  k
  tmp = ft_statfun_depsamplesT(cfg,dat,design(:,permutations(k,:)));
  Ctmp = vismot_cluster6D(cfg2,squareform(tmp.stat));
  pos(k,1) = sum(Ctmp(1).val);
  neg(k,1) = sum(Ctmp(find([Ctmp.thr]<0,1,'first')).val);
end

