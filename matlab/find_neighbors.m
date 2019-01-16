function [output, seed] = find_neighbors(seed, sourcemodel)
%assume seed and sourcemodel in same units
% determine resolution
xpos = sourcemodel.pos(:,1);
resolution = mode(diff(unique(xpos))); % in same units as sourcemodel

if size(seed,2)==1
    % assume seed are indices, covert to x,y,z pos
    seed = sourcemodel.pos(seed,:);
end
[a, b] = size(seed);
output=zeros(a,b,7);
output(:,:,1) = seed;
output(:,:,2) = seed + repmat([-1 0 0], [a,1])*resolution;
output(:,:,3) = seed + repmat([1 0 0], [a,1])*resolution;
output(:,:,4) = seed + repmat([0 -1 0], [a,1])*resolution;
output(:,:,5) = seed + repmat([0 1 0], [a,1])*resolution;
output(:,:,6) = seed + repmat([0 0 -1], [a,1])*resolution;
output(:,:,7) = seed + repmat([0 0 1], [a,1])*resolution;
end