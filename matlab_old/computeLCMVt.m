function [source] = computeLCMVt(source, varargin)

inside  = source.inside;
ninside = length(inside);
npos    = size(source.pos,1);
ntim    = length(source.time);

if length(varargin)==1,
  %do a two-pass single trial analysis correcting for the mean
  doublepass = 1;
  sel        = [1 1];
else
  doublepass = 0;
  sel        = [1 2];
end

filt = cat(1, source.avg.filter{source.inside});
for m = sel
  cnd   = checkdata(varargin{m}, 'datatype', 'timelock');
  nrpt  = size(cnd.trial,1);
  nchan = size(cnd.trial,2);

  dat   = zeros(ninside, ntim);
  datsq = zeros(ninside, ntim); 
  for k = 1:nrpt
    fprintf('analysing single trial data for condition %d, %d/%d\n', m, k, nrpt);
    tmpdat = reshape(cnd.trial(k,:,:),[nchan ntim]);
    tmpdat(~isfinite(tmpdat)) = 0;
    tmpdat = filt*tmpdat;
    if doublepass && exist('dat1', 'var'),
      tmpdat = tmpdat./repmat(sum(dat1.*dof1(ones(size(dat1,1),1),:),2)./sum(dof1), [1 size(dat1,2)])-1;
    end
    dat   = dat   + tmpdat;
    datsq = datsq + tmpdat.^2;
  end  
  dof    = cnd.dof(1,:);
  vardat = (datsq - (dat.^2)./repmat(dof,[size(dat,1) 1]))./(repmat(dof,[size(dat,1) 1])-1);
  dat    = dat./repmat(dof,[size(dat,1) 1]);
  clear datsq;

  if (length(varargin)>1 && m==1) || doublepass,
    dat1 = dat; vardat1 = vardat; dof1 = dof;
  end
end

if length(varargin)>1,
  T = (dat1-dat)./sqrt(vardat1./repmat(dof1,[size(dat1,1) 1])+vardat./repmat(dof,[size(dat,1) 1]));
else
  T = dat./sqrt(vardat./repmat(dof,[size(dat,1) 1]));
end
source.tstat = zeros(npos,ntim);
source.tstat(inside,:) = T;
