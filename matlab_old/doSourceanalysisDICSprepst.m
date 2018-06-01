function [stat,fwhm] = doSourceanalysisDICSprepst(subject, frequency)

%perform T-test on post stimulus versus pre stimulus interval
%collapsed across the four conditions

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
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpst{k}.cumtapcnt);
  ntrl(1,k+4) = length(freqpre{k}.cumtapcnt);
end

oklabel = gradlabel(ok==4);
for k = 1:4
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', frequency);
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

freq{1}   = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
freq{2}   = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
clear freqpst freqpre;
freq      = selectdata(freq{:}, 'param', 'fourierspctrm');

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
%5:8 are the pre stim intervals for conditions 1:4

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
krn  = compute_kernel(sd, 'truncate', 1e-5);
sd   = rmfield(sd, 'fwhm');
sd   = rmfield(sd, 'leadfield');
sd   = rmfield(sd, 'rough');
sd   = rmfield(sd, 'q');

%keep track of this
outside = sd.outside;
inside  = sd.inside;
dim     = pos2dim3d(newgrid.pos);
nchunks = 4;
div     = floor(linspace(0,length(inside),nchunks+1));

cfg.method      = 'pcc';
cfg.keepmom     = 'yes';

%then divide the dipole grid into chunks, to be glued together afterwards
for kk = 1:nchunks  
  cfg.grid.inside = inside((div(kk)+1):div(kk+1));
  %project the single trials through the filters
  tmpsource       = sourceanalysis(cfg, freq);
  tmpsource       = source2sparse(tmpsource);
  tmpsd{kk}       = checkdata(tmpsource, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos', 'inside', 'index');
  clear tmpsource;
end
sd = selectdata(tmpsd{:}, 'param', 'pow');
try,
  sd = rmfield(sd, 'leadfield');
catch
end
clear tmpsd;

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.parameter = 'pow';
cfg2.statistic = 'indepsamplesT';
cfg2.ivar      = 1;

cfg3           = cfg2;
cfg3.statistic = 'yuent';

%congruency T-test
%activation condition 1 vs 3 response left
%activation condition 4 vs 2 response right
conditions = [1 3;4 2];
%for mm = 1:2
    
  %indx1 = btrl(conditions(mm,1)):etrl(conditions(mm,1));
  %indx2 = btrl(conditions(mm,2)):etrl(conditions(mm,2));
  %indx  = [indx1 indx2];
  indx1 = btrl(1):etrl(4);
  indx2 = btrl(5):etrl(8);
  indx  = [indx1 indx2];
  tmpsd = selectdata(sd, 'rpt', indx);
  tmpsd.pow = standardise(log10(tmpsd.pow), 1);
  n1    = numel(indx1);
  n2    = numel(indx2);
  cfg2.design      = [ones(1,n1) 2*ones(1,n2)]; %congruency regressor
  cfg3.design      = cfg2.design;

  tmpstat         = sourcestatistics(cfg2, tmpsd);
  tmpstat.inside  = inside;
  tmpstat.outside = outside;
  tmpstat.pos     = grid.pos;
  dum             = zeros(size(tmpstat.pos,1),1);
  dum(inside)     = tmpstat.stat;
  tmpstat.stat    = dum;
  tmpstat         = rmfield(tmpstat,'mask');
  tmpstat         = rmfield(tmpstat,'prob');
  try, tmpstat         = rmfield(tmpstat,'dim'); end
  tmpstat.stat2   = smooth_vol(tmpstat.stat,krn,dim,inside);
  tmpstat.cfg     = [];
  tmpstaty        = sourcestatistics(cfg3, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.staty   = dum;
  tmpstat.stat2y  = smooth_vol(tmpstat.staty,krn,dim,inside); 
  stat = tmpstat;  
%  if mm==1,
%    stat13 = tmpstat;
%  elseif mm==2,
%    stat42 = tmpstat;
%  end
%end
