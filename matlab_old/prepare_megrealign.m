%first get head motion data

%load volumeconductor
volname = 
vol     = read_vol(volname);


[inda,indb] = match_str({'nasX';'nasY';'nasZ';'lpaX';'lpaY';'lpaZ';'rpaX';'rpaY';'rpaZ'},hmdata.label);
for k = 1:length(hmdata.trial)
  if k==1,
    pnt    = hmdata.trial{k}(indb,:)';
    pnttrl = hmdata.trial{k}(indb,1)';
  else
    pnt    = [pnt; hmdata.trial{k}(indb,:)'];
    pnttrl = [pnttrl; hmdata.trial{k}(indb,1)'];
  end
end

if strcmp(class(pnt), 'single'),
  pnt = double(pnt);
  pnttrl = double(pnttrl);
end
pnt = mean(pnt,1);

T      = headcoordinates(pnt(1:3),pnt(4:6),pnt(7:9));
T(1:3,4) = T(1:3,4)*100; %go to m
tmpvol = transform_vol(inv(T),vol);

cfg         = [];
cfg.vol     = tmpvol;
cfg.grad    = convert_units(hmdata.grad, 'cm');
cfg.channel = hmdata.label;
[tmpvol, grad] = prepare_headmodel(cfg); %singleshell model position specific parameters

cfg.spheremesh  = 642;
cfg.inwardshift = 1;
cfg  = checkconfig(cfg, 'createsubcfg', {'grid'});
grid = prepare_dipole_grid(cfg, tmpvol, grad);
pos  = grid.pos;

allpos = zeros([size(pos) length(hmdata.trial)]);
for k = 1:length(hmdata.trial)
  fprintf('computing dipole positions per trial %d/%d\n',k,length(hmdata.trial));
  tmpT = headcoordinates(pnttrl(k,1:3),pnttrl(k,4:6),pnttrl(k,7:9));
  tmpT(1:3,4) = tmpT(1:3,4)*100;
  tmpvol = transform_vol(inv(tmpT),vol);

  cfg.vol        = tmpvol;
  %[tmpvol, grad] = prepare_headmodel(cfg);
  tmpgrid        = prepare_dipole_grid(cfg, tmpvol, grad);
  allpos(:,:,k)  = tmpgrid.pos;
end

dpos = allpos-repmat(pos,[1 1 length(hmdata.trial)]);
