function [x] = orthogonalise(input,stflag)

if nargin==1,
  stflag = 0;
end

input = standardise(input,1);
x     = input(:,1);

for i = 2:size(input,2)
  y   = input(:,i);
  y   = y-x*((x'*x)\(x'*y));
  if any(y) && stflag, y = standardise(y, 1); end
  if any(y), x = [x y]; end
end