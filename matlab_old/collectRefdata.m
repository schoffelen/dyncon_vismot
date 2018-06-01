function [refdata] = collectRefdata(subject,trl,condition,dftflag)

if nargin==1,
  trl       = [];
  condition = [1:5]; 
  dftflag   = 0;
elseif nargin==2,
  condition = [1:5];
  dftflag   = 0;
elseif nargin==3,
  dftflag   = 0;
end

ncondition = length(condition);
trlnew     = zeros(0,3);
for k = 1:ncondition
  condk = condition(k);
  [dum, tmp{k}] = doPreprocessing(subject,condk,0,1,dftflag);
  if ~isempty(trl),
    trlorig       = findcfg(tmp{k}.cfg, 'trl');
    [int,ia,ib]   = intersect(trl, trlorig, 'rows');
    tmp{k}.trial  = tmp{k}.trial(ib);
    tmp{k}.time   = tmp{k}.time(ib);
    tmp{k}.cfg.trl = trlorig(ib,:);
    trlnew = [trlnew; trl(ia,:)];
  end
end

if length(tmp)>1,
  refdata = appenddata([],tmp{:});
else
  refdata = tmp{1};
end
