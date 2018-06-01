function [stat] = doSourceanalysisLCMV(subject, flag)

if nargin==1,
  flag = 0; %default = unaligned data
end

%-----------------------------------------------------------------
%compute spatial filters based on average single trial covariances
%of all conditions concatenated
cd([subject.pathname,'data/']);
if flag==1,
  load([subject.name,'data_aligned']);
else
  load([subject.name,'data']);
end
data = appenddata([], data1, data2, data3, data4, data5);
warning off; data = struct2double(data); warning on;

%timelockanalysis
cfg              = [];
cfg.vartrllength = 2;
cfg.latency      = [-0.2 0.8];
cfg.covariance   = 'yes';
cfg.blc          = 'yes';
cfg.covariancewindow = [-0.2 0.8];
cfg.blcwindow    = [-0.2 0];
cfg.channel      = 'MEG';
tlck             = timelockanalysis(cfg, data);
clear data;

%load and prune grid
grad  = tlck.grad;
[a,b] = match_str(tlck.label, grad.label);
load([subject.pathname,'grid/',subject.name,'grid6mm']);
for k = 1:length(grid.inside)
  indx  = grid.inside(k);
  tmplf = grid.leadfield{indx};
  grid.leadfield{indx} = tmplf(b,:);
end

%load vol and get it into dewar space
vol       = read_vol([subject.pathname,'vol/',subject.name,'vol.mat']);
if strcmp(subject.datafile(1),'h'),
  %get vol in dewar space
  load([subject.pathname,'dewar2head_avg/',subject.name,'dewar2head_avg.mat'])
  vol       = transform_vol(inv(M), vol);
end

%get grad in correct units
tlck.grad = convert_units(tlck.grad, 'cm');
data1.grad = convert_units(data1.grad, 'cm');
data2.grad = convert_units(data2.grad, 'cm');
data3.grad = convert_units(data3.grad, 'cm');
data4.grad = convert_units(data4.grad, 'cm');

%compute spatial filters and noise
cfg2 = [];
cfg2.method     = 'lcmv';
cfg2.grid       = grid;
cfg2.vol        = vol;
cfg2.fixedori   = 'yes';
cfg2.keepfilter = 'yes';
cfg2.keepmom    = 'yes';
cfg2.projectnoise = 'yes';
cfg2.feedback   = 'none';
cfg2.lambda     = '10%';
source          = sourceanalysis(cfg2,tlck);

filt   = source.avg.filter;
noise  = source.avg.noise;
allmom = source.avg.mom; 
clear source;

%---------------------------------------
%get single conditions in tlck-structure
cfg.keeptrials = 'yes';
cfg.channel    = tlck.label;

warning off;data1 = struct2double(data1);warning on;
data1 = timelockanalysis(cfg, data1);
warning off;data3 = struct2double(data3);warning on;
data3 = timelockanalysis(cfg, data3);
warning off;data2 = struct2double(data2);warning on;
data2 = timelockanalysis(cfg, data2);
warning off;data4 = struct2double(data4);warning on;
data4 = timelockanalysis(cfg, data4);

if length(data1.time)~=length(tlck.time),
  sel = ismember(round(tlck.time.*tlck.fsample),round(data1.time.*tlck.fsample));
  siz = size(data1.trial);
  tmp = zeros([siz(1:2) length(tlck.time)]);
  tmp(:,:,sel) = data1.trial;
  data1.trial  = tmp;
  data1.time   = tlck.time
end
if length(data2.time)~=length(tlck.time),
  sel = ismember(round(tlck.time.*tlck.fsample),round(data2.time.*tlck.fsample));
  siz = size(data2.trial);
  tmp = zeros([siz(1:2) length(tlck.time)]);
  tmp(:,:,sel) = data2.trial;
  data2.trial  = tmp;
  data2.time   = tlck.time
end
if length(data3.time)~=length(tlck.time),
  sel = ismember(round(tlck.time.*tlck.fsample),round(data3.time.*tlck.fsample));
  siz = size(data3.trial);
  tmp = zeros([siz(1:2) length(tlck.time)]);
  tmp(:,:,sel) = data3.trial;
  data3.trial  = tmp;
  data3.time   = tlck.time
end
if length(data4.time)~=length(tlck.time),
  sel = ismember(round(tlck.time.*tlck.fsample),round(data4.time.*tlck.fsample));
  siz = size(data4.trial);
  tmp = zeros([siz(1:2) length(tlck.time)]);
  tmp(:,:,sel) = data4.trial;
  data4.trial  = tmp;
  data4.time   = tlck.time
end

%get single trial time courses (experimental)
for k = 1:size(data1.trial,1)
  data1.cov(k,:,:) = tlck.cov;
end
for k = 1:size(data3.trial,1)
  data3.cov(k,:,:) = tlck.cov;
end
for k = 1:size(data2.trial,1)
  data2.cov(k,:,:) = tlck.cov;
end
for k = 1:size(data4.trial,1)
  data4.cov(k,:,:) = tlck.cov;
end

cfg2.rawtrial     = 'yes';
cfg2.keepfilter   = 'no';
cfg2.projectnoise = 'no';
cfg2.keepmom      = 'yes';

cfgs                  = [];
cfgs.method           = 'montecarlo';
cfgs.statistic        = 'indepsamplesT';
cfgs.numrandomization = 0;
cfgs.parameter        = 'mom';

%-------------------------------------------------
%divide grid in chuncks, do statistics and combine
div  = fix(linspace(0,length(grid.inside),8));
tmpx = zeros(0,length(tlck.time));
tmpy = zeros(0,length(tlck.time));
for k = 1:(length(div)-1)
  tmpgrid         = grid;
  tmpgrid.filter  = cell(1,size(grid.pos,1));
  tmpgrid.inside  = grid.inside((div(k)+1):div(k+1));
  tmpgrid.outside = setdiff(1:size(grid.pos,1), tmpgrid.inside);
  tmpgrid.filter(tmpgrid.inside) = filt(tmpgrid.inside);
  cfg2.grid       = tmpgrid;
  
  tmps1           = sourceanalysis(cfg2, data1);
  tmps1           = source2sparse(tmps1);
  tmps1.dimord    = 'pos';
  tmps1           = checkdata(tmps1, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos_time');
  tmps1           = rmfield(tmps1, 'dim');
  tmps1.dim       = size(tmps1.mom);
  df1             = tmps1.df;

  tmps3           = sourceanalysis(cfg2, data3);
  tmps3           = source2sparse(tmps3);
  tmps3.dimord    = 'pos';
  tmps3           = checkdata(tmps3, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos_time');
  tmps3           = rmfield(tmps3, 'dim');
  tmps3.dim       = size(tmps3.mom);
  df3             = tmps3.df;

  tmps{1}         = tmps1;
  tmps{2}         = tmps3;
  tmps            = selectdata(tmps{:}, 'param', 'mom');
  clear tmps1 tmps3

  krn = gausswin(12,3)'; krn = krn./sum(krn);
  sel = nearest(tlck.time,0);
  for kk = 1:size(tmps.mom,1)
    tmpmom           = squeeze(tmps.mom(kk,:,:));
    tmpmom           = blc(tmpmom, [1 sel]);
    tmps.mom(kk,:,:) = convn(tmpmom, krn, 'same');
  end

  cfgs.design = [ones(1,df1) ones(1,df3)*2];
  tmpstat     = sourcestatistics(cfgs, tmps);

  tmpx        = [tmpx; reshape(tmpstat.stat, [length(tmps.inside) length(tmps.time)])];
  clear tmps tmpstat
  
  tmps2           = sourceanalysis(cfg2, data2);
  tmps2           = source2sparse(tmps2);
  tmps2.dimord    = 'pos';
  tmps2           = checkdata(tmps2, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos_time');
  tmps2           = rmfield(tmps2, 'dim');
  tmps2.dim       = size(tmps2.mom);
  df2             = tmps2.df;

  tmps4           = sourceanalysis(cfg2, data4);
  tmps4           = source2sparse(tmps4);
  tmps4.dimord    = 'pos';
  tmps4           = checkdata(tmps4, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos_time');
  tmps4           = rmfield(tmps4, 'dim');
  tmps4.dim       = size(tmps4.mom);
  df4             = tmps4.df;
  
  tmps{1}         = tmps4;
  tmps{2}         = tmps2;
  tmps            = selectdata(tmps{:}, 'param', 'mom');
  clear tmps2 tmps4

  krn = gausswin(12,3)'; krn = krn./sum(krn);
  for kk = 1:size(tmps.mom,1)
    tmps.mom(kk,:,:) = convn(squeeze(tmps.mom(kk,:,:)), krn, 'same');
  end

  cfgs.design = [ones(1,df4) ones(1,df2)*2];
  tmpstat     = sourcestatistics(cfgs, tmps);

  tmpy        = [tmpy; reshape(tmpstat.stat, [length(tmps.inside) length(tmps.time)])];
  clear tmps tmpstat
end
stat = grid;
stat = rmfield(stat,'leadfield');
stat = rmfield(stat,'dim');
stat.dimord = 'pos_time';
stat.stat13   = zeros(size(stat.pos,1),length(tlck.time));
stat.stat13(stat.inside,:) = tmpx; clear tmpx;
stat.stat42 = zeros(size(stat.pos,1),length(tlck.time));
stat.stat42(stat.inside,:) = tmpy; clear tmpy;
stat.time   = tlck.time;

%compute fwhm for later use
addpath /home/jan/projects/ccc/3D
grid.avg.filter = filt;
grid        = estimate_fwhm_tetra(grid);
stat.fwhm   = grid.fwhm;

%store some additional stuff
stat.noise  = noise;
stat.mom    = zeros(size(stat.stat13));
stat.mom(stat.inside,:) = cat(1,allmom{stat.inside});
