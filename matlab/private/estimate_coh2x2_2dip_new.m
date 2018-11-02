function [output] = estimate_coh2x2_2dip_new(sourcemodel, freq, varargin)

% [COH] = ESTIMATE_COH2x2_2DIP(SOURCEMODEL, FREQ)
%
% estimate the coherence between all dipole-pairs.
% using a pairwise dipole spatial filter

ft_hastoolbox('cellfunction', 1);

lambda    = ft_getopt(varargin, 'lambda',    []);
refindx   = ft_getopt(varargin, 'refindx',   []);
memory    = ft_getopt(varargin, 'memory',    'high');
targetori = ft_getopt(varargin, 'targetori', []);
outputflags = ft_getopt(varargin, 'outputflags',   [1 1 1 1]);

% determine the output requested
assert(isequal(size(outputflags),[1 4]));
getcoh  = outputflags(1);
getcoh0 = outputflags(2); % fitted estimate of 'null'-coherence
getcohc = outputflags(3); % corrected coherence estimate, based on Wens' trick
getcoh1 = outputflags(4); % coherence obtained with a single dipole in the model

% compute cross-spectra
fprintf('computing cross-spectrum for data ...\n');
tmp = freq;
tmp = ft_checkdata(tmp, 'cmbrepresentation', 'fullfast');
C   = tmp.crsspctrm;
clear tmp;

% setting some variables
inside  = sourcemodel.inside; if islogical(inside), inside = find(inside); end
ninside = numel(inside);
nchan   = numel(freq.label);

% reduce to two column leadfields, if needed
nori = size(sourcemodel.leadfield{inside(1)},2);
if nori==3
  fprintf('rotating leadfields into their 2-dimensional basis\n');
  fprintf('and norm normalizing leadfields with the norm of the leadfield\n');
  lf = zeros(nchan, ninside*2);
  for k = 1:ninside
    tmplf   = sourcemodel.leadfield{inside(k)};
    [u,s,v] = svd(tmplf, 'econ');
    tmplf   = tmplf*v(:,1:2);
    lf(:,(k-1)*2+(1:2)) = tmplf./norm(tmplf);
  end
elseif nori==2
  lf = zeros(nchan, ninside*2);
  for k = 1:ninside
    tmplf  = sourcemodel.leadfield{inside(k)};
    lf(:,(k-1)*2+(1:2)) = tmplf./norm(tmplf);
  end
elseif nori==1
  lf = zeros(nchan, ninside);
  for k = 1:ninside
    tmplf   = sourcemodel.leadfield{inside(k)};
    lf(:,k) = tmplf./norm(tmplf);
  end
else
  error('number of orientations unsupported');
end

% % scale the values in the csd matrix; this shouldn't affect the results,
% % but most likely leads to less extreme exponents (numerical issues)
% s1  = svd(real(C));
% C   = C./s1(1);

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(Cf)/size(Cf,1);
else
  lambda = 0;
end

% get a temporary copy for real part, and regularize
rC    = real(C);
rCreg = rC + eye(size(C,1))*lambda;

fprintf('computing some preliminary matrices\n');
lfC     = lf'/rCreg; % lf'*inv(rCreg);
switch memory
  case 'high'
    lfClf   = lfC*lf;    % lf'*inv(rCreg)*lf;
    lfCClfC = lfC*C*lfC';
  case 'low'
end

if nori>1
  fprintf('computing the orientation for the scalar filter\n');
  fprintf('and projecting the leadfields\n');
  newlf = zeros(size(lf,1), size(lf,2)/2);
  ori = zeros(size(newlf,2),2);
  for k = 1:ninside
    indx = (k-1)*2+(1:2);
    switch memory
      case 'high'
        pow  = lfClf(indx,indx)\lfCClfC(indx,indx)/lfClf(indx,indx);
        [u,s,v]    = svd(real(pow));
      case 'low'
        lfClf   = lfC(indx,:)*lf(:,indx);
        lfCClfC = lfC(indx,:)*C*lfC(indx,:)';
        pow     = lfClf\lfCClfC/lfClf;
        [u,s,v] = svd(real(pow));
    end
    newlf(:,k) = normc(lf(:,indx)*u(:,1));
    ori(k,:)   = u(:,1);
  end
  lf = newlf;
elseif nori==1
  ori = ones(size(lf,2),1);
end

fprintf('recomputing some preliminary matrices\n');
lfC     = lf'/rCreg; % lf'*inv(rCreg);
switch memory
  case 'high'
    lfClf   = lfC*lf;    % lf'*inv(rCreg)*lf;
    lfCClfC = lfC*C*lfC';
    lfClfC  = lfC*lfC';
    %lflf    = lf'*lf;
    sumlfClf = sum(lfC.*lf',2); %denominator for single dipole formulation
  case 'low'
    % don't recompute
    lfCC = lfC*C;
    sumlfClf   = sum(lfC.*lf',2);
    sumlfCClfC = sum(lfCC.*conj(lfC),2);
    sumlfClfC  = sum(lfC.*lfC,2);
end

if isempty(refindx)
  refindx = 1:ninside;%100;
end
coh  = zeros(ninside,numel(refindx));
if getcoh0, coh0 = coh; end
if getcohc, cohc = coh; end
if getcoh1, coh1 = coh; end

cnt = 0;
a   = zeros(ninside,1);
r   = zeros(ninside,1);

for k = refindx(:)'
  cnt = cnt+1;
  if mod(k,100)==0, fprintf('computing coherence etc. for voxel %d/%d\n',k,ninside); end
  
  % precompute the diagonals
  lfClf_diag   = diag(lfClf);
  lfCClfC_diag = diag(lfCClfC);
  lfClfC_diag  = diag(lfClfC);
  
  % coherence computation
  switch memory
    case 'high'
      % use pre-computed diagonal matrices for speed up
      denom = inv2x2(convertsquareto2x2(lfClf, k, lfClf_diag));
      denom(:,:,k) = 1./lfClf(k,k);
      numer = convertsquareto2x2(lfCClfC, k, lfCClfC_diag);
      lfClfCtmp = convertsquareto2x2(lfClfC, k, lfClfC_diag);
    case 'low'
      
      % this is the version with the inverse of C sandwiched between lf
      lfClf = zeros(2,2,size(lf,2));
      lfClf(1,2,:) = lfC*lf(:,k);
      lfClf(2,1,:) = conj(lfClf(1,2,:));
      lfClf(1,1,:) = sumlfClf;%sum(lfC.*lf',2);
      lfClf(2,2,:) = lfC(k,:)*lf(:,k);
      
      lfCClfC = zeros(2,2,size(lf,2));
      lfCClfC(1,2,:) = lfCC*lfC(k,:)';
      lfCClfC(2,1,:) = conj(lfCClfC(1,2,:));
      lfCClfC(1,1,:) = sumlfCClfC;%sum(lfCC.*conj(lfC),2);
      lfCClfC(2,2,:) = lfCC(k,:)*lfC(k,:)';
      
      lfClfC = zeros(2,2,size(lf,2));
      lfClfC(1,2,:) = lfC*lfC(k,:)';
      lfClfC(2,1,:) = conj(lfClfC(1,2,:));
      lfClfC(1,1,:) = sumlfClfC;%sum(lfC.*lfC,2);
      lfClfC(2,2,:) = lfC(k,:)*lfC(k,:)';
      
      lfClfCtmp     = lfClfC; % to be used outside the switch
      
      
      denom = inv2x2(lfClf);
      numer = lfCClfC;
  end
  
  if getcohc || getcoh1
    % add (lf'*C*lf)^-1 to denom, as if it were a single dipole scan
    denom(3,1,:) = 1./sumlfClf;
    denom(4,2,:) = denom(3,1,k);
  end
  
  sel   = [1:k-1 k+1:size(lf,2)];
  csd   = sandwichMx2(denom, numer);
  if getcoh0
    csd0 = sandwich2x2(denom(1:2,1:2,:), lfClfCtmp); %filter correlation
    % estimate 'scaling' parameter for null-csd
    [a(cnt),r(cnt)] = fitslope(permute(real(csd0(1,1,sel)),[3 1 2]),permute(real(csd(1,1,sel)),[3 1 2]));
    %r(cnt) = 0;
    %a(cnt) = median(permute(csd0(1,2,sel),[3 1 2])./permute(real(csd(1,2,sel)),[3 1 2]));
  end
  
  pp = sqrt(abs(csd(1,1,:)).*abs(csd(2,2,:)));
  coh(:,cnt)  = csd(1,2,:)./pp;
  coh(k,cnt)  = 1;
  
  if getcoh0
    coh0(sel,cnt) = (csd0(1,2,sel).*a(cnt))./pp(sel);
    coh0(k,  cnt) = 1;
  end
  
  if getcohc
    %csd = sandwich2x2(denom([1 4],:,:), numer);
    pp = sqrt(abs(csd(1,1,:)).*abs(csd(4,4,:)));
    cohc(:,cnt)   = csd(1,4,:)./pp;
  end
  
  if getcoh1
    %csd = sandwich2x2(denom([3 4],:,:), numer);
    pp = sqrt(abs(csd(3,3,:)).*abs(csd(4,4,:)));
    coh1(:,cnt) = csd(3,4,:)./pp;
    cohc(k,cnt) = 0;
  end
end

fprintf('creating output structure\n');
if getcoh,  output.coh   = single(coh);  clear coh;  end
if getcoh0, output.coh0  = single(coh0); clear coh0; end
if getcohc, output.cohc  = single(cohc); clear cohc; end
if getcoh1, output.coh1  = single(coh1); clear coh1; end
output.ori   = ori;

if isempty(targetori)
  targetori = repmat([1 0 0],ninside,1);
end

if isfield(sourcemodel, 'v')
  signvec = zeros(ninside,1);
  for k = 1:ninside
    ik = inside(k);
    this = sourcemodel.v{ik}*output.ori(k,:)';
    signvec(k) = sign(this(:)'*targetori(k,:)'+eps);
  end
  
  if getcoh,  output.coh  = output.coh.*(signvec*signvec(refindx)');  end
  if getcoh0, output.coh0 = output.coh0.*(signvec*signvec(refindx)'); end
  if getcohc, output.cohc = output.cohc.*(signvec*signvec(refindx)'); end
  if getcoh1, output.coh1 = output.coh1.*(signvec*signvec(refindx)'); end
  output.ori  = output.ori.*repmat(signvec,[1 size(output.ori,2)]);
end

output.a = a;
output.r = r;

%output.lf = lf;

function [a,r,m,n] = fitslope(x,y)

xorig = x;

x = abs(x);
y = abs(y);

%[t1, t2] = percthreshold(x, 0.5, 1, 0);
t1  = prctile(x, 95);
sel = x>t1;

x  = x(sel);
mx = mean(x);

y  = y(sel);
my = mean(y);

%w = spdiags(x(:), 0, n, n);
x = x-mx;
y = y-my;
%a = (x'*w*y)./(x'*w*x);
xsq = x.^2;
a = (xsq'*y)./(xsq'*x); % this is the same as the weighted regression above, with weights the values of x
r = sum((y-a*x).^2)./(y'*y);

if nargout>2
  m = a*xorig;
  n = sum(sel);
end