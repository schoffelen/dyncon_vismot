function vismot_connectivitybrowser(grid, source, parameter, parameter2, anasc, cohsc)

% HCP_CONNECTIVITYBROWSER allows to browse a 6D-connectome
% 
% Use as
%  
% hcp_connectivitybrowser(grid, source, parameter)
%
% where grid is the specification of a 3D grid that specified the
% sourcemodel for the connectome, and source is the structure containing
% the connectome (in the sparse representation, i.e. only the inside
% voxels), and parameter is a string denoting the fieldname in the source
% structure denoting the connectivity metric

h = figure;

if nargin<5
  anasc = [];
end
if nargin<6
  cohsc = [];
end

coh = zeros(prod(grid.dim), 'single')+nan;
coh2 = zeros(prod(grid.dim), 'single')+nan;

if ~isa(source.(parameter), 'single')
  source.(parameter) = single(source.(parameter));
end
if ~isa(source.(parameter2), 'single')
  source.(parameter2) = single(source.(parameter2));
end
coh(grid.inside, grid.inside) = source.(parameter);
coh = reshape(coh, [grid.dim grid.dim]);
coh2(grid.inside, grid.inside) = source.(parameter2);
coh2 = reshape(coh2, [grid.dim grid.dim]);


ana = zeros(grid.dim);
ana(grid.inside) = .5;

siz = size(coh);
x1i = round(siz(1)/2);
y1i = round(siz(2)/2);
z1i = round(siz(3)/2);
x2i = round(siz(4)/2);
y2i = round(siz(5)/2);
z2i = round(siz(6)/2);

x = 1:siz(1);
y = 1:siz(2);
z = 1:siz(3);

% ensure same color scaling for all figures
if isempty(anasc)
  c1min = nanmin(coh(:));
  c1max = nanmax(coh(:));
else
  c1min = anasc(1);
  c1max = anasc(2);
end

% ensure same color scaling for all figures
if isempty(cohsc)
  c2min = nanmin(coh2(:));
  c2max = nanmax(coh2(:));
else
  c2min = cohsc(1);
  c2max = cohsc(2);
end
if c2max==c2min
  c2max = c2max+10000*eps(1);
end
if isempty(c2min) && isempty(c2max)
  c2min = 0;
  c2max = 1;
end

val1 = 0;
val2 = 0;
while ishandle(h)
  clf
  h1 = subplot( 'position', [0.01 0.15 0.48 0.82]); %axis([1 siz(2) 1 siz(3)]); % subplot(2,4,1);
  h2 = subplot( 'position', [0.51 0.15 0.48 0.82]); %axis([1 siz(1) 1 siz(3)]); % subplot(2,4,2);
  h1t = subplot('position', [0.01 0.02 0.48 0.12]);
  h2t = subplot('position', [0.51 0.02 0.48 0.12]);
  
  vol1 = shiftdim(coh(:,:,:,x2i,y2i,z2i));
  vol2 = shiftdim(coh2(x1i,y1i,z1i,:,:,:));
  
  [map1,siz1,ndiv1] = vol2map(vol1);
  [map2,siz2,ndiv2] = vol2map(vol2);
    
  subplot(h1); 
  hmap1 = imagesc(map1);axis xy;axis off;caxis([c1min c1max]); 
  [nx,ny] = ind2sub(ndiv1,z2i);
  ystar   = siz1(2)*(ny-1)+y2i;
  xstar   = siz1(1)*(nx-1)+x2i;
  hold on;plot(xstar,ystar,'w*','linewidth',2); set(gca, 'tag', 'left');
  colormap jet;
  
  subplot(h2); 
  hmap2 = imagesc(map2);axis xy;axis off;caxis([c2min c2max]); 
  [nx,ny] = ind2sub(ndiv2,z1i);
  ystar   = siz2(2)*(ny-1)+y1i;
  xstar   = siz2(1)*(nx-1)+x1i;
  hold on;plot(xstar,ystar,'w*','linewidth',2); set(gca, 'tag', 'right');
  colormap jet
  
  ind1 = sub2ind(grid.dim,x1i,y1i,z1i);
  ind2 = sub2ind(grid.dim,x2i,y2i,z2i);
  
  
  p1   = grid.pos(ind1,:);
  p2   = grid.pos(ind2,:);
  subplot(h1t);text(0,.5,sprintf('ind1=%6.0f\npos1=[%3.1f %3.1f %3.1f],\nval1=%f', ind1,p1(1),p1(2),p1(3), val1));axis off
  subplot(h2t);text(0,.5,sprintf('ind2=%6.0f\npos2=[%3.1f %3.1f %3.1f],\nval2=%f', ind2,p2(1),p2(2),p2(3), val2));axis off

  try, [d1, d2, key] = ginput(1); catch, return; end
  if key==113  % q
    return
  end
  % update the view to a new position
  t  = get(gca, 'tag');
  
  switch t
  case 'left'
    % update the indices for the volume to be displayed on the right
    [x1i,y1i,z1i] = two2three(d1,d2,size(map1),siz1);
    val1 = vol1(x1i,y1i,z1i);
  case 'right'
    % update the indices for the volume to be displayed on the right
    [x2i,y2i,z2i] = two2three(d1,d2,size(map2),siz2);
    val2 = vol2(x2i,y2i,z2i);
  end
  x1i = round(x1i);
  y1i = round(y1i);
  z1i = round(z1i);
  x2i = round(x2i);
  y2i = round(y2i);
  z2i = round(z2i);
  if x1i<1, x1i=1; end
  if y1i<1, y1i=1; end
  if z1i<1, z1i=1; end
  if x2i<1, x2i=1; end
  if y2i<1, y2i=1; end
  if z2i<1, z2i=1; end
  
  
    if x1i>siz(1), x1i=siz(1); end
    if y1i>siz(2), y1i=siz(2); end
    if z1i>siz(3), z1i=siz(3); end
    if x2i>siz(1), x2i=siz(1); end
    if y2i>siz(2), y2i=siz(2); end
    if z2i>siz(3), z2i=siz(3); end
end

function [map,siz,ndiv] = vol2map(vol)

siz    = size(vol);
ndiv   = [ceil(sqrt(siz(3))) floor(sqrt(siz(3)))];
map    = zeros(siz(2)*ndiv(2),siz(1)*ndiv(1));
for k=1:siz(3)
  [nx,ny] = ind2sub(ndiv,k);
  map(siz(2)*(ny-1)+1:siz(2)*ny,siz(1)*(nx-1)+1:siz(1)*nx) = vol(:,:,k)';
end

function [x,y,z] = two2three(d1,d2,siz2d,siz3d)

% siz2d(1) is an integer multiple of siz3d(2)
% siz2d(2) is an integer multiple of siz3d(1)
% slices start in the lower left corner

ndiv1 = siz2d(2)./siz3d(1);
ndiv2 = siz2d(1)./siz3d(2);

d1 = round(d1);
d2 = round(d2);

mapx = zeros(siz3d(2)*ndiv2,siz3d(1)*ndiv1);
mapy = mapx;
mapz = mapx;
for k=1:siz3d(3)
  [nx,ny] = ind2sub([ndiv1 ndiv2],k);
  mapx(siz3d(2)*(ny-1)+(1:siz3d(2)),siz3d(1)*(nx-1)+(1:siz3d(1))) = repmat(1:siz3d(1),[siz3d(2) 1]);
  mapy(siz3d(2)*(ny-1)+(1:siz3d(2)),siz3d(1)*(nx-1)+(1:siz3d(1))) = repmat(1:siz3d(2),[siz3d(1) 1])';
  mapz(siz3d(2)*(ny-1)+(1:siz3d(2)),siz3d(1)*(nx-1)+(1:siz3d(1))) = k;
end
x = mapx(d2,d1);
y = mapy(d2,d1);
z = mapz(d2,d1);
