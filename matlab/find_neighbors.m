function [output, seed] = find_neighbors(seed, sourcemodel)
%assume seed and sourcemodel in same units
% determine resolution
xpos = sourcemodel.pos(:,1);
resolution = mode(diff(unique(xpos))); % in same units as sourcemodel

if size(seed,2)==1
    % assume seed are indices, covert to x,y,z pos
    seed = sourcemodel.pos(seed,:);
end

output=[];
for k=1:size(seed,1)
    output(end+1,:) = seed(k,:);
    output(end+1,:) = seed(k,:) + [-1 0 0]*resolution;
    output(end+1,:) = seed(k,:) + [1 0 0]*resolution;
    output(end+1,:) = seed(k,:) + [0 -1 0]*resolution;
    output(end+1,:) = seed(k,:) + [0 1 0]*resolution;
    output(end+1,:) = seed(k,:) + [0 0 -1]*resolution;
    output(end+1,:) = seed(k,:) + [0 0 1]*resolution;
end