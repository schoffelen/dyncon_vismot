function [stat, design] = doSourceanalysisDICSroicoh(subject, frequency)

%Try to do the something similar as doSourceanalysisDICSpst2

fieldtripdefs

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

cd([subject.pathname,'freq/']);
if frequency<37,
  load([subject.name,'mtmfft004']);
else
  load([subject.name,'mtmfft012']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
rtall = [];
blockpst = [];
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpst{k}.cumtapcnt);
  trl   = findcfg(freqpst{k}.cfg, 'trl');
  rt    = trl(:,4);
  rtall = [rtall; rt];
  freqpst{k}.rt = rt./fsample;
  %blockpst = [blockpst trl2blockregressor(subject, trl)];
end
block = blockpst;

oklabel = gradlabel(ok==4);
for k = 1:4
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', frequency);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
end

freq    = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
freq.rt = rtall;
clear freqpst

%load grid and prune leadfields
[a,b] = match_str(freq.label, gradlabel);
if flag<2,
  load([subject.pathname,'grid/',subject.name,'grid6mm.mat']);
elseif flag==2,
  load([subject.pathname,'grid/',subject.name,'grid6mmNew.mat']);
end
eval('newgrid = grid;');
for k = 1:length(newgrid.inside)
  indx  = newgrid.inside(k);
  tmplf = newgrid.leadfield{indx};
  newgrid.leadfield{indx} = tmplf(b,:);
end

vol       = read_vol([subject.pathname,'vol/',subject.name,'vol.mat']);
if strcmp(subject.datafile(1),'h'),
  %get vol in dewar space
  load([subject.pathname,'dewar2head_avg/',subject.name,'dewar2head_avg.mat'])
  vol       = transform_vol(inv(M), vol);
end

%get grad in correct units
freq.grad = convert_units(freq.grad, 'cm');

cfg            = [];
cfg.method     = 'dics';
cfg.frequency  = frequency;
cfg.grid       = newgrid;
cfg.vol        = vol;
cfg.fixedori   = 'yes';
cfg.keepfilter = 'yes';
cfg.realfilter = 'yes';
cfg.keepleadfield = 'yes';
cfg.powmethod  = 'lambda1';
cfg.lambda     = '10%';
cfg.feedback   = 'none';
source         = sourceanalysis(cfg, freq);

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);
%1:4 are the post stim intervals for conditions 1:4

%project the leadfields onto the optimal orientation
cfg.grid.filter = source.avg.filter;
for k = 1:length(source.inside)
  indx = source.inside(k);
  cfg.grid.leadfield{indx} = newgrid.leadfield{indx}*source.avg.ori{indx};
end

%compute FWHM of spatial filters
source.leadfield = cfg.grid.leadfield;
addpath /home/jan/projects/ccc/3D
sd   = estimate_fwhm_tetra(source);
fwhm = sd.fwhm;
krn  = compute_kernel(sd, 'truncate', 2e-5);
sd   = rmfield(sd, 'fwhm');
sd   = rmfield(sd, 'leadfield');
sd   = rmfield(sd, 'rough');
sd   = rmfield(sd, 'q');

%keep track of this
outside = sd.outside;
inside  = sd.inside;
numel(inside)
dim     = pos2dim3d(newgrid.pos);
nchunks = 4;
div     = floor(linspace(0,length(inside),nchunks+1));

cfg.method      = 'pcc';
cfg.keepmom     = 'yes';
tmpsource       = ft_sourceanalysis(cfg, freq);
keyboard
%then divide the dipole grid into chunks, to be glued together afterwards
for kk = 1:nchunks  
  cfg.grid.inside = inside((div(kk)+1):div(kk+1));
  %project the single trials through the filters
  tmpsource       = ft_sourceanalysis(cfg, freq);
  tmpsource       = source2sparse(tmpsource);
  tmpsd{kk}       = checkdata(tmpsource, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos');
  clear tmpsource;
end
sd = selectdata(tmpsd{:}, 'param', 'pow');
try,
  sd = rmfield(sd, 'leadfield');
catch
end
clear tmpsd;

%load('/home/jan/projects/visuomotor/matlab/alphaROI');
%load('/home/jan/projects/visuomotor/matlab/erdROI');
load(['/analyse/4/Project0030/roi/',subject.name,'roimaskInducedGamma1'],'roi1','roi2');
x1 = find(ismember(find(isfinite(fwhm)), roi1));
x2 = find(ismember(find(isfinite(fwhm)), roi2));
%x1 = find(ismember(find(isfinite(fwhm)), sel1));
%x2 = find(ismember(find(isfinite(fwhm)), sel2));

roi1 = standardise(nanmean(log10(sd.pow(:,x1)),2),1);
roi2 = standardise(nanmean(log10(sd.pow(:,x2)),2),1);

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.parameter = 'pow';
cfg2.statistic = 'glm';
cfg2.glm.statistic = 'T';
cfg2.glm.contrast  = [0 0 1 0];

cfg3              = cfg2;
cfg3.glm.contrast = [0 0 0 1];

cfg4              = cfg2;
cfg4.glm.contrast = [0 0 0 0 0 1 0 0];

cfg5              = cfg2;
cfg5.glm.contrast = [0 0 0 0 0 0 0 1]; 

%glm T contrast for regressor with ROI activity,
%activation condition 1 vs 3 response left congruent vs incongruent
%activation condition 4 vs 2 response right congruent vs incongruent
%activation condition 1 vs 2 visual stimulus left congruent vs incongruent
%activation condition 4 vs 3 visual stimulus right congruent vs incongruent
conditions = [1 2 3 4];
sd.dim = size(sd.pow);
for mm = 1:size(conditions,1)
    
  indx1 = btrl(conditions(mm,1)):etrl(conditions(mm,1));
  indx2 = btrl(conditions(mm,2)):etrl(conditions(mm,2));
  indx3 = btrl(conditions(mm,3)):etrl(conditions(mm,3));
  indx4 = btrl(conditions(mm,4)):etrl(conditions(mm,4));
  indx  = [indx1 indx2 indx3 indx4];
  tmpsd = selectdata(sd, 'rpt', indx);
  tmpsd.pow = standardise(log10(tmpsd.pow), 1);
  tmpsd.dim = size(tmpsd.pow);
  n1    = numel(indx1);
  n2    = numel(indx2);
  n3    = numel(indx3);
  n4    = numel(indx4);
  cfg2.design      = [ones(1,n1) -ones(1,n2) -ones(1,n3)  ones(1,n4)]; %congruency    regressor
  cfg2.design(2,:) = [ones(1,n1) -ones(1,n2)  ones(1,n3) -ones(1,n4)]; %response side regressor
  cfg2.design(3,:) = roi1(indx)';
  cfg2.design(4,:) = roi2(indx)';
  tmp = orthogonalise(cfg2.design([1 2 3],:)',1)'; cfg2.design(3,:) = tmp(3,:);
  tmp = orthogonalise(cfg2.design([1 2 4],:)',1)'; cfg2.design(4,:) = tmp(3,:);
  cfg3.design      = cfg2.design;
  design{mm,1}     = cfg2.design;
  design{mm,2}     = cfg3.design;

  tmpstat         = sourcestatistics(cfg2, tmpsd);
  tmpstat.inside  = inside;
  tmpstat.outside = outside;
  tmpstat.pos     = grid.pos;
  dum             = zeros(size(tmpstat.pos,1),1);
  dum(inside)     = tmpstat.stat;
  tmpstat.stat    = dum;
  tmpstat         = rmfield(tmpstat,'mask');
  tmpstat         = rmfield(tmpstat,'prob');
  tmpstat         = rmfield(tmpstat,'dim');
  tmpstat.stat2   = smooth_vol(tmpstat.stat,krn,dim,inside);
  tmpstat.cfg     = [];
  tmpstaty        = sourcestatistics(cfg3, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statb   = dum;
  tmpstat.stat2b  = smooth_vol(tmpstat.statb,krn,dim,inside); 
  
  %fields stat and stat2 are T-statistic of correlation with roi1
  %fields statb and stat2b are T-statistic of correlation with roi2

  lcli = [ones(1,n1)  zeros(1,n2) -ones(1,n3) zeros(1,n4)];
  lcri = [ones(1,n1)  -ones(1,n2) zeros(1,n3) zeros(1,n4)];
  rcli = [zeros(1,n1) zeros(1,n2) -ones(1,n3)  ones(1,n4)];
  rcri = [zeros(1,n1) -ones(1,n2) zeros(1,n3)  ones(1,n4)];

  lcli = lcli - mean(lcli);
  lcri = lcri - mean(lcri);
  rcli = rcli - mean(rcli);
  rcri = rcri - mean(rcri);

  roi1 = roi1(indx)' - mean(roi1(indx));
  roi2 = roi2(indx)' - mean(roi2(indx));

  ppi1resp = roi1.*rcri;
  ppi1vis  = roi1.*rcli;
  ppi2resp = roi2.*lcli;
  ppi2vis  = roi2.*lcri;

  cfg4.design      = cfg2.design;
  cfg4.design(3,:) = roi1;
  cfg4.design(4,:) = roi2;
  cfg4.design(5,:) = (rcri); %abs(rcri);
  cfg4.design(6,:) = ppi1resp;
  cfg4.design(7,:) = (lcli);
  cfg4.design(8,:) = ppi2resp;
  
  tmpdesign = cfg4.design;
  tmp  = orthogonalise(tmpdesign([1:4 5 6 7 8],:)',1)'; cfg4.design([7 8],:) = tmp(7:8,:);
  tmp  = orthogonalise(tmpdesign([1:4 7 8 5 6],:)',1)'; cfg4.design([5 6],:) = tmp(7:8,:);
  design{mm,3} = cfg4.design;
  
  cfg5.design  = cfg4.design;
  cfg5.design(3,:) = roi1;
  cfg5.design(4,:) = roi2;
  cfg5.design(5,:) = (rcli);
  cfg5.design(6,:) = ppi1vis;
  cfg5.design(7,:) = (lcri);
  cfg5.design(8,:) = ppi2vis;

  tmpdesign = cfg5.design;
  tmp  = orthogonalise(tmpdesign([1:4 5 6 7 8],:)',1)'; cfg5.design([7 8],:) = tmp(7:8,:);
  tmp  = orthogonalise(tmpdesign([1:4 7 8 5 6],:)',1)'; cfg5.design([5 6],:) = tmp(7:8,:);
  design{mm,4} = cfg5.design;
  
  %ppi regressors based on congruency contrast (same response side congruent minus incongruent)
  cfg4.glm.contrast = [0 0 0 0 0 1 0 0];
  tmpstaty        = sourcestatistics(cfg4, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statppi1resp  = dum;
  tmpstat.stat2ppi1resp = smooth_vol(tmpstat.statppi1resp,krn,dim,inside);
  cfg4.glm.contrast = [0 0 0 0 0 0 0 1];
  tmpstaty        = sourcestatistics(cfg4, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statppi2resp  = dum;
  tmpstat.stat2ppi2resp = smooth_vol(tmpstat.statppi2resp,krn,dim,inside);
  cfg4.glm.contrast = [0 0 0 0 0 1 0 -1];
  tmpstaty        = sourcestatistics(cfg4, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statdppiresp  = dum;
  tmpstat.stat2dppiresp = smooth_vol(tmpstat.statdppiresp,krn,dim,inside);
  
  %ppi regressors based on informative visual stimulus contrast (same hemifield processing target stimulus)
  cfg5.glm.contrast = [0 0 0 0 0 1 0 0];
  tmpstaty        = sourcestatistics(cfg5, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statppi1vis  = dum;
  tmpstat.stat2ppi1vis = smooth_vol(tmpstat.statppi1vis,krn,dim,inside);
  cfg5.glm.contrast = [0 0 0 0 0 0 0 1];
  tmpstaty        = sourcestatistics(cfg5, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statppi2vis  = dum;
  tmpstat.stat2ppi2vis = smooth_vol(tmpstat.statppi2vis,krn,dim,inside);
  cfg5.glm.contrast = [0 0 0 0 0 1 0 -1];
  tmpstaty        = sourcestatistics(cfg5, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.statdppivis  = dum;
  tmpstat.stat2dppivis = smooth_vol(tmpstat.statdppivis,krn,dim,inside);
  tmpstat.cfg     = [];
  tmpstat.dim     = pos2dim3d(tmpstat.pos);

  if mm==1,
    stat = tmpstat;
  end
end
