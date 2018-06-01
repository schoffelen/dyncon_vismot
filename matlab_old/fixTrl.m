function [newtrl] = fixTrl(subject, trl)

%get a trl-matrix which is consistent with the representation in 
%subjinfo (according to analyzeRT2). This means that for sessions
%consisting of multiple runs the values in the first two columns
%of trl will be adjusted according to some offset.
%In general the trl matrix can be meaningless due to redefenition
%of resampled trials. Be aware of this!


tmp  = getEvent(subject);

for k = 1:numel(trl)
  begsmp = trl{k}(:,1);
  dsmp   = diff([1000000;begsmp]);
  nblock(1,k) = sum(dsmp<0);
  bblock{1,k} = find(dsmp<0);
end

if any(nblock ~= numel(tmp)),
  error('discrepancy in number of runs');
end

offset = zeros(1,numel(tmp));
for k = 1:numel(tmp)-1
  offset(1,(k+1):end) = offset(1,(k+1):end) + tmp{k}(end).sample + 1000;
end

newtrl = zeros(0,size(trl{1},2));
for k = 1:numel(trl)
  tmptrl = trl{k}; 
  ix     = bblock{k};
  ix     = [ix;size(tmptrl,1)+1];
  for m = 1:numel(ix)-1
    tmptrl(ix(m):(ix(m+1)-1),1:2) = tmptrl(ix(m):(ix(m+1)-1),1:2) + offset(m);
  end
  newtrl = [newtrl; tmptrl];
end
