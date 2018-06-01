%get a grad structure
grad = convert_units(grad,'cm');
pnt  = grad.pnt(1:248,:);

[O,r] = fitsphere(pnt);
dpnt  = pnt-O(ones(1,248),:);


sel = 1:248;
nchan = length(sel);
%go to spherical coordinates
[Phi,Th,R] = cart2sph(dpnt(sel,1),dpnt(sel,2),dpnt(sel,3));
Th         = Th+0.5*pi;
cosTh      = cos(Th)';
Phi        = Phi';
R          = R';

lmax = 8;

%get basis functions
nbasis = (lmax+1).^2-1;
lmindx = zeros(nbasis,2);
endsmp = 0;
for llop = 1:lmax
  begsmp = endsmp+1;
  endsmp = endsmp+(2*llop+1);
  lmindx(begsmp:endsmp,1) = llop;
  lmindx(begsmp:endsmp,2) = -llop:llop;
end

onevecC = ones(1,nchan);
onevecB = ones(nbasis,1);

Y = zeros(nbasis,nchan) + i.*zeros(nbasis,nchan);
for llop = 1:lmax
  l = llop;
  m = [-l:l]';
  N = sqrt([(2*l+1)/(4*pi)].*[factorial(l-m)./factorial(l+m)]);
  mpos = [0:l]';
  P    = legendre(l, cosTh).*[(-1).^mpos(:,onevecC)];
  mneg = -m(1:l);
  tmp  = [(-1).^mneg.*factorial(l-mneg)./factorial(l+mneg)];
  %tmp  = [(-1).^mneg];
  Pneg = tmp(:,onevecC).*P(end:-1:2,:);
  P    = [Pneg;P]; 
  expPhi = exp(i.*(m*Phi));
  Y(lmindx(:,1)==llop,:) = P.*expPhi.*N(:,onevecC);
end

a = transpose(Y./(R(onevecB,:).^(lmindx(:,onevecC)+1)));
b = transpose(Y.*(R(onevecB,:).^(lmindx(:,onevecC)  )));
