%needed variables: hdr, mri, shape
coils = hdr.orig.user_block_data{12}.pnt(6:10,:); %positions in headspace
coils = coils*100;

%%convert units
%shape = convert_units(shape, 'mm');
%
%%extract scalp surface
%scp   = scalp3b(mri.anatomy);
%scp   = warp_apply(mri.transform, scp);
%nscp  = size(scp,1);
%
%%estimate closest scalp points
%for k = 1:size(coils,1)
%  tmp  = scp - coils(k*ones(nscp,1),:);
%  dtmp = sqrt(sum(tmp.^2,2));
%
%end

%why not do something like this:
%wobble around head position (tolerance 1.5 in each direction)
%wobble around rotation (tolerance 20 degrees in each direction)
%get a set of positions so that the leadfields need to be computed only once

for xtrans = -5:5
for ytrans = -5:5
for ztrans = -10:5
for xrot = -10:10
  for yrot = -10:10
    for zrot = -10:10
      x   = xrot*pi/180;
      y   = yrot*pi/180;
      z   = zrot*pi/180;
      rx  = [1 0 0;0 cos(x) -sin(x);0 sin(x) cos(x)];
      ry  = [cos(y) 0 sin(y);0 1 0;-sin(y) 0 cos(y)];
      rz  = [cos(z) -sin(z) 0;sin(z) cos(z) 0;0 0 1];
      R   = rz*ry*rx;
      T   = [xtrans;ytrans;ztrans];
      M   = [R T;0 0 0 1];
      %pos{xtrans+6,ytrans+6,ztrans+11,xrot+11,yrot+11,zrot+11} = (R*coils')'-coils;
      pos{xtrans+6,ytrans+6,ztrans+11,xrot+11,yrot+11,zrot+11} = warp_apply(M, coils)-coils;
    end
  end
end
end
end
end


%----some datasets contain an initial COH-run
%this allows for an estimate of both initial position, and orientation of the coils
%load in header of datafile with COH runs
pnt = hdr.orig.user_block_data{15}.pnt;
ori = hdr.orig.user_block_data{15}.ori;


%----try something else
%compute position per trial
hmdata = headmotiontracking(cfg, data);

%compute transformation matrix per trial
coils = [cfg.coils ones(5,1)];
M     = cell(1,length(hmdata));
for k = 1:length(hmdata)
  tmpM = [hmdata{k}.dip.allpos ones(5,1)]\coils;
  tmpM = tmpM';
  tmpM(4,:) = [0 0 0 1];
  M{k} = tmpM;
end

grad = cell(1,length(hmdata));
for k = 1:length(hmdata)
  grad{k} = transform_sens(M{k}, data.grad);
end

