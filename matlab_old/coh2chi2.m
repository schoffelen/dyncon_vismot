function [output] = coh2chi2(dat,dof)

if size(dat,3) ~= numel(dof)
  error('numel(dof) should be equal to size(dat,3)');
end

%transform dat (assuming dat to be unsquared coherence)
dat = atanh(abs(dat));

%make dof same size as dat (memory wise not smart but convenient for product)
doforig = dof;
dof = zeros(1,1,numel(doforig));
dof(:) = doforig;
dof = repmat(dof, [size(dat(:,:,1)) 1]);

output = 2.*(sum(dof.*dat.^2,3) - (sum(dof.*dat,3).^2)./sum(dof,3));

