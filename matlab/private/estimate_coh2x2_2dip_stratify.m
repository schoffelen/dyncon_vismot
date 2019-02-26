function [output] = estimate_coh2x2_2dip_stratify(sourcemodel, freq, varargin)

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
stratflag = ft_getopt(varargin, 'stratflag', true);

% determine the output requested
assert(isequal(size(outputflags),[1 4]));
getcoh  = outputflags(1);
getcoh0 = outputflags(2); % fitted estimate of 'null'-coherence
getcohc = outputflags(3); % corrected coherence estimate, based on Wens' trick
getcoh1 = outputflags(4); % coherence obtained with a single dipole in the model

% compute cross-spectra
if ~iscell(freq) && numel(freq)~=2
  error('input freq should be a 1x2 cell-array');
end
fprintf('computing cross-spectrum for data ...\n');
[C, n] = freq2C(freq);

if ~stratflag
  C1 = freq2C(freq{1});
  C2 = freq2C(freq{2});
end


% setting some variables
inside  = sourcemodel.inside; if islogical(inside), inside = find(inside); end
ninside = numel(inside);
nchan   = numel(freq{1}.label);

if isempty(refindx)
  refindx = 1:ninside;%100;
end

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

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  regu = ratio * trace(C)/size(C,1);
elseif isempty(lambda)
  regu = 0;
else
  regu = lambda;
end

% get a temporary copy for real part, and regularize
rC    = real(C);
rCreg = rC + eye(size(C,1))*regu;

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

if stratflag>0
  lfC     = lf'/rCreg; % lf'*inv(rCreg);
  sumlfClf = sum(lfC.*lf',2); %denominator for single dipole formulation

  % single dipole model spatial filters
  w = spdiags(1./sumlfClf,0,ninside,ninside)*lfC;
  
  fprintf('computing single trial power\n');
  [pow, n, P] = getpow(freq, w(refindx,:)); % this is done in a subfunction
  pow = log10(pow);
  
  % stratify for power
  powindx = dostratification(pow,n, 10);
  
  
  if stratflag==1
    % these variables will be used below
    powindx1 = sparse(double(       powindx(:,1:n(1))))*P(       1:n(1),1:sum(freq{1}.cumtapcnt)); % indexes the rows in the fourier matrix
    powindx2 = sparse(double(powindx(:,n(1)+(1:n(2)))))*P(n(1)+(1:n(2)),(sum(freq{1}.cumtapcnt)+1):end);
  elseif stratflag==2
    n1 = n(1);
    n2 = n(2);
    
    while n1~=n2
      rr = sum(powindx)./size(powindx,1);
      if n1<n2
        % we may want to equalize the number of trials, by throwing out data
        % from condition 2
        rr(1:n(1)) = nan;
        rr(rr==0) = nan;
        [~,ix] = nanmin(rr);
        
      elseif n2<n1
        % we may want to equalize the number of trials, by throwing out data
        % from condition 1
        rr(n(1)+(1:n(2))) = nan;
        rr(rr==0) = nan;
        [~,ix] = nanmin(rr);
        
      end
      pow(:,ix) = nan;
      powindx   = dostratification(pow,n, 10);
      
      n1 = n(1)-double(sum(sum(powindx(:,1:n(1)))==0));
      n2 = n(2)-double(sum(sum(powindx(:,n(1)+(1:n(2))))==0));
    end
    
    rr   = sum(powindx)./size(powindx,1);
    sel1 = find(rr(1:n(1))~=0);
    sel2 = find(rr(n(1)+(1:n(2)))~=0);
    
    cfg = [];
    cfg.trials = sel1; 
    freq{1} = ft_selectdata(cfg, freq{1});
    cfg.trials = sel2;
    freq{2} = ft_selectdata(cfg, freq{2});
    
    [C,n] = freq2C(freq);
    C1    = freq2C(freq{1});
    C2    = freq2C(freq{2});
        
    % it is difficult to give a quantitative estimate of lambda, therefore also
    % support relative (percentage) measure that can be specified as string (e.g. '10%')
    if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
      ratio = sscanf(lambda, '%f%%');
      ratio = ratio/100;
      regu = ratio * trace(C)/size(C,1);
    elseif isempty(lambda)
      regu = 0;
    else
      regu = lambda;
    end
    
    % get a temporary copy for real part, and regularize
    rC    = real(C);
    rCreg = rC + eye(size(C,1))*regu;

  end
end

% recompute some stuff, now the orientation has been fixed, and possibly
% the number of trials stratified (stratflag==2)
fprintf('recomputing some preliminary matrices\n');
lfC     = lf'/rCreg; % lf'*inv(rCreg);
switch memory
  case 'high'
    lfClf   = lfC*lf;    % lf'*inv(rCreg)*lf;
    lfCClfC = lfC*C*lfC';
    lfClfC  = lfC*lfC';
    %lflf    = lf'*lf;
    sumlfClf = sum(lfC.*lf',2); %denominator for single dipole formulation
    
    if stratflag~=1
      lfCC1lfC = lfC*C1*lfC';
      lfCC2lfC = lfC*C2*lfC';
    end
    
    % precompute the diagonals
    lfClf_diag   = diag(lfClf);
    lfCClfC_diag = diag(lfCClfC);
    lfClfC_diag  = diag(lfClfC);
    
    if stratflag~=1
      lfCC1lfC_diag = diag(lfCC1lfC);
      lfCC2lfC_diag = diag(lfCC2lfC);
    end
    
  case 'low'
    lfCC       = lfC*C;
    sumlfClf   = sum(lfC.*lf',2);
    sumlfCClfC = sum(lfCC.*conj(lfC),2);
    sumlfClfC  = sum(lfC.*lfC,2);
    if ~stratflag
      lfCC1 = lfC*C1;
      lfCC2 = lfC*C2;
      sumlfCC1lfC = sum(lfCC1.*conj(lfC),2);
      sumlfCC2lfC = sum(lfCC2.*conj(lfC),2);
    end
end

coh  = zeros(ninside,numel(refindx));
coh_1 = coh;
coh_2 = coh;
if getcoh0, coh0 = coh; end
if getcohc, cohc = coh; end
if getcoh1, coh1 = coh; end

cnt = 0;
a   = zeros(ninside,1);
r   = zeros(ninside,1);

for k = refindx(:)'
  cnt = cnt+1;
  if mod(k,100)==0, fprintf('computing coherence etc. for voxel %d/%d\n',k,ninside); end
    
  % coherence computation
  switch memory
    case 'high'
      
      if stratflag==1, error('this does not deal with type1 stratification yet'); end
     
      % use pre-computed diagonal matrices for speed up
      denom = inv2x2(convertsquareto2x2(lfClf, k, lfClf_diag));
      denom(:,:,k) = 1./lfClf(k,k);
      numer = convertsquareto2x2(lfCClfC, k, lfCClfC_diag);
      lfClfCtmp = convertsquareto2x2(lfClfC, k, lfClfC_diag);
    
      numer1 = convertsquareto2x2(lfCC1lfC, k, lfCC1lfC_diag);
      numer2 = convertsquareto2x2(lfCC2lfC, k, lfCC2lfC_diag);
      
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
      
      % create the condition specific numerators based on the stratified
      % power. Note that the spatial filters are implicitly computed using the
      % unstratified data
      if stratflag==1
        F1 = transpose(freq{1}.fourierspctrm(powindx1(cnt,:)>0,:));
        F2 = transpose(freq{2}.fourierspctrm(powindx2(cnt,:)>0,:));
        n1 = sum(powindx(cnt,1:n(1)));
        n2 = sum(powindx(cnt,n(1)+(1:n(2))));
        
        lfCC1 = lfC*((F1*F1')./n1);
        lfCC2 = lfC*((F2*F2')./n2);
      end
      lfCC1lfC = zeros(2,2,size(lf,2));
      lfCC1lfC(1,2,:) = lfCC1*lfC(k,:)';
      lfCC1lfC(2,1,:) = conj(lfCC1lfC(1,2,:));
      lfCC1lfC(2,2,:) = lfCC1(k,:)*lfC(k,:)';
      
      lfCC2lfC = zeros(2,2,size(lf,2));
      lfCC2lfC(1,2,:) = lfCC2*lfC(k,:)';
      lfCC2lfC(2,1,:) = conj(lfCC2lfC(1,2,:));
      lfCC2lfC(2,2,:) = lfCC2(k,:)*lfC(k,:)';
      
      if stratflag
        lfCC1lfC(1,1,:) = sum(lfCC1.*conj(lfC),2);
        lfCC2lfC(1,1,:) = sum(lfCC2.*conj(lfC),2);
      else
        lfCC1lfC(1,1,:) = sumlfCC1lfC;
        lfCC2lfC(1,1,:) = sumlfCC2lfC;
      end
      
      numer1          = lfCC1lfC;
      numer2          = lfCC2lfC;
      
  end
  
  if getcohc || getcoh1
    % add (lf'*C*lf)^-1 to denom, as if it were a single dipole scan
    denom(3,1,:) = 1./sumlfClf;
    denom(4,2,:) = denom(3,1,k);
  end
  
  sel   = [1:k-1 k+1:size(lf,2)];
  csd   = sandwichMx2(denom, numer);
  
  % this is using a 'glm' type heuristic to estimate the 'null'-coherence,
  % i.e. the stuff that one would expect based in the filter correlations
  % in combination with a single parameter estimate for the noise leakage
  if getcoh0
    csd0 = sandwich2x2(denom(1:2,1:2,:), lfClfCtmp); %filter correlation
    % estimate 'scaling' parameter for null-csd
    [a(cnt),r(cnt)] = fitslope(permute(real(csd0(1,1,sel)),[3 1 2]),permute(real(csd(1,1,sel)),[3 1 2]));
    %r(cnt) = 0;
    %a(cnt) = median(permute(csd0(1,2,sel),[3 1 2])./permute(real(csd(1,2,sel)),[3 1 2]));
  end
  
  % this is computing the normal coherence
  pp = sqrt(abs(csd(1,1,:)).*abs(csd(2,2,:)));
  coh(:,cnt)  = csd(1,2,:)./pp;
  coh(k,cnt)  = 1;
  
  csd_1 = sandwichMx2(denom, numer1); pp_1 = sqrt(abs(csd_1(1,1,:)).*abs(csd_1(2,2,:)));
  csd_2 = sandwichMx2(denom, numer2); pp_2 = sqrt(abs(csd_2(1,1,:)).*abs(csd_2(2,2,:)));
  
  coh_1(:,cnt) = csd_1(1,2,:)./pp_1;
  coh_2(:,cnt) = csd_2(1,2,:)./pp_2;
  
  if getcoh0
    coh0(sel,cnt) = (csd0(1,2,sel).*a(cnt))./pp(sel);
    coh0(k,  cnt) = 1;
  end

  % this is computing the coherence as per Wens' geometric trick
  if getcohc
    %csd = sandwich2x2(denom([1 4],:,:), numer);
    pp = sqrt(abs(csd(1,1,:)).*abs(csd(4,4,:)));
    cohc(:,cnt)   = csd(1,4,:)./pp;
  end
  
  % this is computing the coherence according to a single-dipole source
  % model
  if getcoh1
    %csd = sandwich2x2(denom([3 4],:,:), numer);
    pp = sqrt(abs(csd(3,3,:)).*abs(csd(4,4,:)));
    coh1(:,cnt) = csd(3,4,:)./pp;
    cohc(k,cnt) = 0;
  end
end

fprintf('creating output structure\n');
if getcoh
  output.coh   = single(coh);  clear coh;
  output.coh_1 = single(coh_1); clear coh_1;
  output.coh_2 = single(coh_2); clear coh_2;
end
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
t1  = prctile(x, 90);
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


function powindx = dostratification(pow, n, numbin)

cfgx = [];
cfgx.numbin = numbin;
cfgx.equalbinavg = 'no';
powindx = pow>0;
fprintf('performing the stratification\n');
for k = 1:size(pow,1)
  out = stratify(cfgx,pow(k,1:n(1)),pow(k,n(1)+(1:n(2))));
  powindx(k,:) = [isfinite(out{1}) isfinite(out{2})];
end
  
function [pow, n, P] = getpow(freq, w)

n(1) = numel(freq{1}.cumtapcnt);
n(2) = numel(freq{2}.cumtapcnt);

sumtapcnt = cumsum([0;freq{1}.cumtapcnt;freq{2}.cumtapcnt]);

pow = zeros(size(w,1),sum(n));
for k = 1:sum(n)
  if k<=n(1)
    freqindx = 1;
    indx     = (sumtapcnt(k)+1):(sumtapcnt(k+1));
    indx2    = (sumtapcnt(k)+1):(sumtapcnt(k+1));
  else
    freqindx = 2;
    indx     = ((sumtapcnt(k)+1):(sumtapcnt(k+1))) - sum(freq{1}.cumtapcnt);
    indx2    = ((sumtapcnt(k)+1):(sumtapcnt(k+1)));
  end
  
  pow(:,k)  = mean(abs(w*transpose(freq{freqindx}.fourierspctrm(indx,:))).^2,2);
  P(k,indx2) = 1;
end
P = sparse(P);

function [C, n] = freq2C(freq)

if ~iscell(freq)
  freq = {freq};
end

for k = 1:numel(freq)
  tmp = freq{k};
  tmp = ft_checkdata(tmp, 'cmbrepresentation', 'fullfast');
  C(:,:,k) = tmp.crsspctrm;
  n(1,k)   = numel(freq{k}.cumtapcnt);
  clear tmp;
end

for k = 1:size(C,3)
  if k==1
    tmp = zeros(size(C,1),size(C,1));
  end
  tmp = tmp + C(:,:,k).*n(k);
end
C = tmp./sum(n);

function [varargout] = stratify(cfg, varargin)

% FT_STRATIFY tries to reduce the variance in a specific feature in the data
% that is not related to an effect in two or multiple conditions, but where
% that feature may confound the analysis. Stratification is implemented by
% randomly removing elements from the data, making the distribution of the
% data equal on that feature.
%
% Use as
%   [output]          = ft_stratify(cfg, input1, input2, ...), or
%   [output, binaxis] = ft_stratify(cfg, input1, input2, ...)
%
% For the histogram and the split method, each input is a Nchan X Nobs
% matrix. The output is a cell-array with in each cell the same data as in
% the corresponding input, except that the observations that should be
% removed are marked with a NaN.
%
% For the equatespike method, each input is a Ntrials X 1 cell-array. Each
% trial should contain the spike firing moments (i.e. a logical Nchans X
% Nsamples matrix). The output is a cell-array with in each cell the same
% data as in the corresponding input, except that spike numbers have been
% equated in each trial and channel.
%
% The configuration should contain
%   cfg.method      = 'histogram'
%                     'splithilo'
%                     'splitlohi'
%                     'splitlolo'
%                     'splithihi'
%                     'equatespike'
%
% The following options apply only to histogram and split methods.
%   cfg.equalbinavg = 'yes'
%   cfg.numbin      = 10
%   cfg.numiter     = 2000
%
% The following options apply only to the equatespike method.
%   cfg.pairtrials  = 'spikesort', 'linkage' or 'no' (default = 'spikesort')
%   cfg.channel     = 'all' or list with indices ( default = 'all')
%
% See also FT_FREQSTATISTICS, FT_TIMELOCKSTATISTICS, FT_SOURCESTATISTICS

% input1 and input2 are the to be stratified with respect to each other
% dimensionality of input1 (2) = chan x rpt. If nchan>1, do a "double"
% stratification

% set the defaults
cfg.equalbinavg  = ft_getopt(cfg, 'equalbinavg', 'yes');
cfg.numbin       = ft_getopt(cfg, 'numbin', 10);
cfg.numiter      = ft_getopt(cfg, 'numiter', 2000);
cfg.binedges     = ft_getopt(cfg, 'binedges', []);

% the input data is a cell-array containing matrices for each condition
input = varargin;

if size(input{1},1)~=size(input{2},1)
  ft_error('the number of channels should be the same');
end
if size(input{1},1)~=1 && strcmp(cfg.equalbinavg, 'yes')
  ft_error('combining equalising bin-averages and simultaneous stratification of two channels is impossible');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan = size(input{1},1);
ncond = length(input);
% if nchan~=2, ft_error('number of trials ~= 2, do not know how to stratify'); end

if ~isempty(cfg.binedges)
  fprintf('using the bin edges for the histogram as defined in the cfg\n');
  if iscell(cfg.binedges)
    cfg.numbin = cellfun(@numel, cfg.binedges)-1;
  else
    cfg.numbin = numel(cfg.binedges)-1;
  end
end

%------prepare some stuff
if numel(cfg.numbin) ~= nchan
  cfg.numbin = repmat(cfg.numbin(1), [1 nchan]);
end
for j = 1:nchan
  tmp  = [];
  for cndlop = 1:ncond
    tmp = cat(2, tmp, input{cndlop}(j,:));
  end % concatenate input-data
  
  % create a 'common binspace'
  [ndum,x] = hist(tmp, cfg.numbin(j));
  dx    = diff(x);
  if ~isempty(cfg.binedges)
    if iscell(cfg.binedges)
      xc = cfg.binedges{j};
    else
      xc = cfg.binedges;
    end
  else
    xc    = [-inf x(1:end-1)+0.5*dx(1) inf];
  end
  
  %make histograms and compute the 'target histogram', which
  %will equalize the conditions while throwing away as few
  %repetitions as possible
  tmp = [];
  for cndlop = 1:ncond
    [n{cndlop}(j,1:numel(xc)), b{cndlop}(j,:)] = histc(input{cndlop}(j,:), xc);
    tmp = [tmp; n{cndlop}(j,:)];
  end
  binaxis(j,1:(cfg.numbin(j)+1)) = xc;
end

%------index the trials
%------create n-d histo
linearhisto = zeros(ncond, prod(cfg.numbin));
for cndlop = 1:ncond
  tmpb = zeros(1, size(b{cndlop},2));
  for j = 1:nchan
    if j == 1
      tmpb = tmpb + (b{cndlop}(j,:)).*prod(cfg.numbin(1:(j-1)));
    else
      tmpb = tmpb + (b{cndlop}(j,:)-1).*prod(cfg.numbin(1:(j-1)));
    end
  end
  b{cndlop}             = tmpb;
  for binlop = 1:size(linearhisto,2)
    linearhisto(cndlop,binlop) = sum(tmpb==binlop);
  end
end

%------find intersection
numok  = min(linearhisto,[],1);
nummax = max(linearhisto,[],1);
for j = 1:length(numok)
  minind{j} = find(linearhisto(:,j)==numok(j));
end

%------do stratification the easy or the hard way
if strcmp(cfg.equalbinavg, 'yes')
  %---this is the hard way
  if nchan>1, ft_error('the option equalbinavg only works for one channel input at the moment'); end
  
  %---first allocate some memory
  for cndlop = 1:ncond
    sel{cndlop}     = zeros(1,size(b{cndlop},2));
  end
  
  for binlop = 1:length(numok)
    if numok(binlop)>0
      tmpmatmin    = nan(ncond,nummax(binlop));
      tmpmatmax    = nan(ncond,nummax(binlop));
      tmpmatminind = nan(ncond,nummax(binlop));
      tmpmatmaxind = nan(ncond,nummax(binlop));
      for cndlop = 1:ncond
        tmpsel          = find(b{cndlop}==binlop);
        tmpdat          = input{cndlop}(tmpsel);
        [tmpsrt,tmpind] = sort(tmpdat);
        tmpmatmin(   cndlop,1:linearhisto(cndlop,binlop)        ) = tmpsrt;
        tmpmatmax(   cndlop,end-linearhisto(cndlop,binlop)+1:end) = tmpsrt;
        tmpmatminind(cndlop,1:linearhisto(cndlop,binlop)        ) = tmpind;
        tmpmatmaxind(cndlop,end-linearhisto(cndlop,binlop)+1:end) = tmpind;
      end
      refavg = nanmean(tmpmatmin,2);
      refmin = mean(tmpmatmin(:,1:numok(binlop)        ),2);
      refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
      
      %---determine the amount of trials in this bin, so that for all conditions
      %it lies between the lowest possible and highest possible realisation
      ok     = 0; lowok = 0; hiok  = 0; cnt   = 0;
      offset = zeros(1,ncond);
      while ok==0
        if numok(binlop)>0
          [tmpref,refind] = min(refavg(minind{binlop}));
          if any(refmin - tmpref > 0)
            numok(binlop)                          = numok(binlop) - 1;
            offset(minind{binlop}(refind))         = offset(minind{binlop}(refind)) + 1; %correction term
            tmpmatmin(minind{binlop}(refind),:)    = [   tmpmatmin(minind{binlop}(refind),2:end) nan];
            %tmpmatminind(minind{binlop}(refind),:) = [tmpmatminind(minind{binlop}(refind),2:end) nan];
            warning off;
            refavg = nanmean(tmpmatmin,2);
            refmin = mean(tmpmatmin(:,1:numok(binlop)),2);
            refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
            warning on;
          else
            lowok = 1;
          end
          [tmpref,refind] = min(refavg(minind{binlop}));
          if any(refmax - tmpref < 0)
            numok(binlop)                          = numok(binlop) - 1;
            tmpmatmax(minind{binlop}(refind),:)    = [nan    tmpmatmax(minind{binlop}(refind),1:end-1)];
            %tmpmatmaxind(minind{binlop}(refind),:) = [nan tmpmatmaxind(minind{binlop}(refind),1:end-1)];
            warning off;
            refavg = nanmean(tmpmatmax,2);
            refmin = mean(tmpmatmin(:,1:numok(binlop)),2);
            refmax = mean(tmpmatmax(:,end-numok(binlop)+1:end),2);
            warning on;
          else
            hiok = 1;
          end
        end
        if lowok==1 && hiok==1, ok = 1; end
        if cnt==100,            ok = 1; end
        cnt = cnt+1;
      end
      
      %---this is now the average that should be approximated
      [tmpref,refind] = min(refavg(minind{binlop}));
      
      if numok(binlop)>0
        for cndlop = 1:ncond
          pnt     = tmpmatmin(cndlop, 1:linearhisto(cndlop,binlop)) - tmpref;
          nrow    = length(pnt)-numok(binlop)+1;
          pntmat  = repmat(pnt,[nrow 1]);
          % get a good initial guess
          cpnt    = conv2(1,ones(1,numok(binlop))./numok(binlop),pnt,'same');
          [~,indc]= min(abs(cpnt));
          indvec  = [indc-floor(numok(binlop)/2):indc+ceil(numok(binlop)/2)-2];
          if length(indvec)<=1,       indvec = [indc indc];                    end
          if indvec(1)<1,             indvec = indvec-indvec(1)+1;             end
          if indvec(end)>length(cpnt),  indvec = indvec-indvec(end)+length(cpnt);  end
          tmpmat  = zeros(nrow,length(pnt));
          tmpmat(:, indvec                       ) = 1;
          if length(unique(indvec))>1 || size(pntmat,2)>nrow
            tmpmat(:, setdiff(1:length(pnt),indvec)) = eye(nrow);
          else
            tmpmat = eye(nrow);
          end
          %tmpmat  = [ones(nrow,numok(binlop)-1) eye(nrow)];
          %tmpmat  = tmpmat(:,randperm(size(tmpmat,2)));
          if cndlop~=minind{binlop}(refind)
            m      = nan(1,100);
            for rndlop = 1:100
              if rndlop<=12 || sum(diff(m(rndlop-11:rndlop-1))==0)<10
                dif = abs(sum(pntmat.*tmpmat,2)./numok(binlop));
                [m(rndlop),ind] = min(dif);
                tmpvec           = tmpmat(ind,:);
                tmpmat           = repmat(tmpvec,[nrow 1]);
                indone           = find(tmpmat(1,:));
                tmpsel           = randperm(length(indone));
                tmpsel           = indone(tmpsel(1));
                tmpmat(:,[tmpsel find(tmpmat(1,:)==0)]) = eye(nrow);
              else
                %do nothing and go on
                break
              end
            end
          else
            tmpvec = [ones(1,numok(binlop)) zeros(1,nrow-1)];
          end
          tmpsel = find(b{cndlop}==binlop);
          sel{cndlop}(tmpsel(tmpmatminind(cndlop,tmpvec))) = 1;
        end
      end
    end
  end
  
else
  %------stratify the easy way
  for cndlop = 1:ncond
    sel{cndlop} = zeros(1,size(b{cndlop},2));
    tmphisto    = linearhisto(cndlop,:);
    for binlop = 1:size(tmphisto,2)
      tmpsel = find(b{cndlop}==binlop);
      tmpsel = tmpsel(randperm(length(tmpsel)));
      tmpsel = tmpsel(1:numok(binlop));
      sel{cndlop}(tmpsel) = 1;
    end
  end
end

%------create output
output = input;
for cndlop = 1:ncond
  output{cndlop}(:,find(sel{cndlop}==0))=nan;
end
varargout{1} = output;
varargout{2} = binaxis;

