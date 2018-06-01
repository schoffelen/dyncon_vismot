function [err, J] = errorfunction(params, dat, sens, vol, constr, reducerank, normalize, normalizeparam)

if strcmp(class(params),'struct'),
  param = params.param;
  dat   = params.dat;
  sens  = params.sens;
  vol   = params.vol;
  constr= params.constr;
  reducerank = params.reducerank;
  normalize = params.normalize;
  normalizeparam = params.normalizeparam;
elseif nargin>1,
  param = params;
else
  error('function requires input as a structure, or as a list of arguments');
end

%param refers to the xyz-coordinates of the origin, and the rotation angles 
%of the fixed array of coils with respect to the sensor coordinate frame

%extract rotation
x  = param(end-2);
y  = param(end-1);
z  = param(end);
rx = [1 0 0;0 cos(x) -sin(x);0 sin(x) cos(x)];
ry = [cos(y) 0 sin(y);0 1 0;-sin(y) 0 cos(y)];
rz = [cos(z) -sin(z) 0;sin(z) cos(z) 0;0 0 1];
R  = rz*ry*rx;

%extract position
n   = size(constr.dpos,1);
pos = repmat(param(1:3), [n 1]) + (R*constr.dpos')';

% construct the leadfield matrix for all dipoles
lf = compute_leadfield(pos, sens, vol, 'reducerank', reducerank, ...
      'normalize', normalize, 'normalizeparam', normalizeparam);
if isfield(constr, 'ori'),
  lf = lf * ori;
end

% compute the optimal dipole moment and the model error
% ordinary least squares, this is the same as MLE with weight=eye(nchans,nchans)
mom = pinv(lf)*dat;
dif = dat - lf*mom;
err = sum(dif(:).^2) / sum(dat(:).^2);

if nargout>1,
  incr = [mean(abs(constr.dpos(:)))*0.01*ones(1,3) 0.005 0.005 0.005];
  J = jacobianJM(@errorfunction, params, incr); 
end
