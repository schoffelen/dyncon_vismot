function pow = fourier2pow(filter, fourierspctrm, cumtapcnt)
% this function creates power values in source space from the combination
% of the fourierspectrum and the spatial filter. The input cumtapcnt
% contains information about amount of tapers per trial.

nrpt = numel(cumtapcnt);
ntap_all = sum(cumtapcnt);

ix = [];
iz = [];
for k=1:nrpt
    ntap = cumtapcnt(k);
    ix = [ix; repmat(k, [ntap 1])];
    iz = [iz; ones(ntap,1)./ntap];
end
iy = 1:ntap_all;

P  = sparse(iy,ix,iz,ntap_all,nrpt);

pow = (abs(filter*transpose(fourierspctrm)).^2)*P;