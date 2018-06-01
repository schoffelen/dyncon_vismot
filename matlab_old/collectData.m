function [data,trlnew] = collectData(subject, trl, condition, dftflag)

%[data] = collectData(subject, trl)
%collect all data across conditions for specified subject. 
%if second input trl is specified collect a subset of trials

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
  [tmp{k}, dum] = doPreprocessing(subject,condk,1,0,dftflag);
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
  data = appenddata([],tmp{:});
else 
  data = tmp{1};
end
