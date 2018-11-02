function pow = fourier2pow(w,fourier,tapvec)

siz = [size(w) 1];
nchan = siz(2);
ndip  = siz(3);
nrpt  = numel(tapvec);
ntap  = size(fourier,1);

if nchan~=size(fourier,2)
  error('unexpected dimensions');
end

% create the projection matrix from tapers to trials
ix = (1:ntap)';
iy = zeros(ntap,1);
iz = zeros(ntap,1);
indx = 0;
for k = 1:numel(tapvec)
  indx = max(indx)+(1:tapvec(k));
  iy(indx) = k;
  iz(indx) = 1./tapvec(k);
end
P = sparse(ix,iy,iz);

w = permute(w, [3 2 1]);
pow = zeros(size(w,1), nrpt);
fourier = transpose(fourier);
for k = 1:size(w,3)
  % sum all components of w
  pow = pow + (abs(w(:,:,k)*fourier).^2)*P;
end
