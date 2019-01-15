function [output, seed] = find_neighbors(seed, sourcemodel)
%assume seed and sourcemodel in same units
% determine resolution
xpos = sourcemodel.pos(:,1);
ref = xpos(1);
nonzero_pos = abs(ref-xpos)>0;
res = min(abs(xpos(nonzero_pos))); % in same units as sourcemodel


if size(seed,2)==1
    % assume seed are indices, covert to x,y,z pos
    seed = sourcemodel.pos(seed,:);
end

output=[];
for k=1:size(seed,1)
    output(end+1,:) = seed(k,:);
    output(end+1,:) = seed(k,:) + [-1 0 0]*res;
    output(end+1,:) = seed(k,:) + [1 0 0]*res;
    output(end+1,:) = seed(k,:) + [0 -1 0]*res;
    output(end+1,:) = seed(k,:) + [0 1 0]*res;
    output(end+1,:) = seed(k,:) + [0 0 -1]*res;
    output(end+1,:) = seed(k,:) + [0 0 1]*res;
end