function [new_refindx, n_neighbors] = revise_neighbors(refindx, insidepos, resolution)

    new_refindx = [];
    n_neighbors = [];
    for l = 1:7:numel(refindx)
      seed_orig = insidepos(refindx(l),:);
      neighb = insidepos(refindx(l+1:l+6),:);
      distance = sqrt(sum((neighb-seed_orig).^2,2));
      nondirect_neighb = find(round(10*distance)~=round(10*resolution));
      tmp = refindx(l:l+6);
      tmp(nondirect_neighb+1) = []; % exclude non-direct neighbors
      tmp = unique(tmp); % exclude multiple copies of same neighbor.
      new_refindx = [new_refindx; tmp];
      n_neighbors = [n_neighbors; numel(tmp)-1];   % determine the amount of valid neighbors per roi
    end
