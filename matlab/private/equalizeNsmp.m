function [indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2)

%this function takes two vectors of numbers (not necessarily of different lengths)
%and returns two vectors containing indices to the elements in the input and two
%vectors containing the number of samples such that a subset from smp1 and smp2
%match in terms of number of elements and number of samples. the smallest number
%of elements and samples wins

smp1 = smp1(:);
smp2 = smp2(:);

n1      = numel(smp1);
n2      = numel(smp2);
[n, in] = min([n1 n2]);

%convert to shortest and longest vector
if in==1,
  smpa = smp1;
  smpb = smp2;
elseif in==2,
  smpa = smp2;
  smpb = smp1;
end

indxa = logical(zeros(size(smpa)));
indxb = logical(zeros(size(smpb)));
allsa = zeros(n,1);
allsb = zeros(n,1);
tmpb  = zeros(0,2);
while any(isfinite(smpa)) && size(tmpb,1)+n<length(smpb),
  [m, im] = min([min(smpa) min(smpb)]);
  if im==1,
    %shortest trial left belongs to shortest sample
    [p, ip]   = min(smpa);
    indxa(ip) = 1;
    allsa(ip) = p;

    [p2, ip2] = min(smpb);
    indxb(ip2) = 1;
    allsb(ip2) = p;

    smpa(ip)  = inf;
    smpb(ip2) = inf;
  elseif im==2,
    %shortest trial left belongs to longer sample
    %we can theoretically do without this trial but 
    %this is not sure yet
    [p2, ip2] = min(smpb);
    tmpb      = [tmpb; ip2 p2];
    smpb(ip2) = inf;
  end 
end

%now if the shortest trial left belongs to the longer sample we have to prune the trial
%in the shorter sample
while any(isfinite(smpa))
  [m, im] = min([min(smpa) min(smpb)]);
  if im==1,
    %shortest trial left belongs to shortest sample
    [p, ip]   = min(smpa);
    indxa(ip) = 1;
    allsa(ip) = p;

    [p2, ip2] = min(smpb);
    indxb(ip2) = 1;
    allsb(ip2) = p;

    smpa(ip)  = inf;
    smpb(ip2) = inf;
  elseif im==2,
    %shortest trial left belongs to longer sample
    [p2, ip2]  = min(smpb);
    indxb(ip2) = 1;
    allsb(ip2) = p2;
    
    [p, ip]    = min(smpa);
    indxa(ip)  = 1;
    allsa(ip)  = p2;

    smpa(ip)  = inf;
    smpb(ip2) = inf;
  end 
end

if in==1,
  indx1 = find(indxa);
  nsmp1  = allsa(indx1);
  indx2 = find(indxb);
  nsmp2  = allsb(indx2);
elseif in==2,
  indx1 = find(indxb);
  nsmp1  = allsb(indx1);
  indx2 = find(indxa);
  nsmp2  = allsa(indx2);
end
