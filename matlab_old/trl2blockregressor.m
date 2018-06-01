function [design, trl] = trl2blockregressor(subject, trl)

cd(subject.pathname);
cd('block');
load([subject.name,'block']);

tmp  = trl(:,1);
dtmp = diff([0;tmp]);
if any(dtmp<0),
  %multiple runs present
  newrun = find(dtmp<0);
  newtrl = trl;
  for k = 1:numel(newrun)
    %add a constant both to block, and trl
    newtrl(newrun(k):end,   1:2) = newtrl(newrun(k):end,   1:2) + 10000000;
    block(2, find(block(1,:)>k)) = block(2, find(block(1,:)>k)) + 10000000;
  end
else
  %only one run present
  newtrl = trl;
end

trl      = newtrl;
trl(:,5) = 0;
tmp      = trl(:,1);
for k = 1:size(block,2)
  sel = find(tmp>block(2,k));
  trl(sel,5) = k;
end

design = zeros(size(block,2)-1, size(trl,1));
for k = 1:size(design,1)
  sel = find(trl(:,5)==k);
  design(k,sel) = 1;
end
