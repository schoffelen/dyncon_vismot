function coh = average_seed_over_neighbors(coh, refindx, n_neighbors, index_orig_seed)

tmp1 = coh.coh;
tmp2 = zeros(size(tmp1,1), numel(n_neighbors));
    index=1;
    for m=1:numel(n_neighbors)
        tmp2(:,m) = nanmean(tmp1(:,index:index+n_neighbors(m)), 2);
        index = (index+n_neighbors(m))+1;
    end
    coh.coh = tmp2;
    % manually set coherence at original refindx to one.
    for m=1:numel(index_orig_seed)
        coh.coh(refindx(index_orig_seed(m)),m) = 1;
    end
end