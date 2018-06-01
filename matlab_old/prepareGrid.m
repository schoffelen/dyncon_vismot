function [grid] = prepareGrid(subject, resolution, grad)

if nargin==1 || isempty(resolution),
  resolution = 6;
end

spm_defaults
defaults.analyze.flip = 0;

%chmt
chmt = strcmp(subject.datafile(1), 'h');

%load mri
cd([subject.pathname,'mri']);
load([subject.name,'mri']);

%load vol
cd([subject.pathname,'vol']);
vol = read_vol([subject.name,'vol.mat']);

if nargin==2 || isempty(grad),
  %load grad
  cd([subject.pathname,'data']);
  load([subject.name,'data'], 'data1');
  grad  = convert_units(data1.grad, 'cm');
else
  grad = convert_units(grad, 'cm');
end
warning off;
grad  = struct2double(grad);
tra     = grad.tra;
balance = grad.balance;
warning on;

if chmt,
  %get dewar2head transformation matrix
  cd([subject.pathname,'dewar2head_avg']);
  load([subject.name,'dewar2head_avg']);
  
  %transform in average dewar space
  vol  = transform_vol(inv(M), vol);
  tmpM = M;
  tmpM(1:3,4)   = tmpM(1:3,4)*10; %mri is in mm
else
  tmpM = eye(4);
  
  %get average grad positions and oris from the headers of the different runs
  tmpname = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname];
  for k = 1:length(subject.runnames)
    tmphdr = read_header([tmpname,subject.runnames{k},'c,rfDC']);
    if k==1,
      gradorig = tmphdr.grad;
      cohorig  = tmphdr.orig.user_block_data{15}.pnt(1:3,:);
      coh      = tmphdr.orig.user_block_data{15}.pnt(1:3,:);
      cohorig  = cohorig + tmphdr.orig.user_block_data{17}.pnt(1:3,:);
      coh      = coh     + tmphdr.orig.user_block_data{17}.pnt(1:3,:);
    else
      coh      = coh + tmphdr.orig.user_block_data{15}.pnt(1:3,:);
      coh      = coh + tmphdr.orig.user_block_data{17}.pnt(1:3,:);
    end
  end
  %these are the average coil-on-head locations
  cohorig = cohorig./2;
  coh     = coh./(2*k);

  %construct some meaningful transformation matrix (exact meaning is irrelevant);
  Morig = headcoordinates(cohorig(3,:),cohorig(1,:),cohorig(2,:));
  M     = headcoordinates(coh(3,:),    coh(1,:),    coh(2,:));
  T     = M*inv(Morig); %to transform grad belonging to first run into average grad
  grad  = transform_sens(T, gradorig);
  grad  = convert_units(grad, 'cm');
  grad.tra     = tra;
  grad.balance = balance;
end

%load template grid
try
  fname = [homejanmatlab,'/mri/templategrid',num2str(resolution),'mm'];        
  load(fname, 'grid');
catch
  error('no templategrid of this resolution exists, you first have to create one');
end

cfg             = [];
cfg.coordinates = 'ctf';
cfg.nonlinear   = 'no';
%cfg.template    = [homejanmatlab,'/spm2/templates/T1.mnc'];
cfg.template    = [homejanmatlab,'/toolboxes/spm2/templates/T1.mnc'];
normalise       = volumenormalise(cfg,mri);
pos             = warp_apply(inv(tmpM)*inv(normalise.cfg.final), grid.pos*10, 'homogenous')/10;
grid.pos        = pos;

newgrad = grad;
%newgrad.tra(1:248,249:end) = 0;
cfg             = [];
cfg.grid        = grid;
cfg.grad        = newgrad;
cfg.vol         = vol;
cfg.channel     = channelselection('MEG',grad.label);
grid            = prepare_leadfield(cfg);

%save grid
cd([subject.pathname,'grid']);
fname = [subject.name,'grid',num2str(round(resolution)),'mm'];        
if exist([fname,'.mat'],'file')
  fname = [fname,'New'];
end
save(fname, 'grid');
