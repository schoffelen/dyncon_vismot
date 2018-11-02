function [y, x_norm] = normc(x)

% NORMC column-wise norm normalisation

x_norm = sqrt(sum(x.^2));
y      = x*diag(1./x_norm);
