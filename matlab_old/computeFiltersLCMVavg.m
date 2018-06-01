function [source] = computeFiltersLCMVavg(subject, flag)

if nargin==1,
  flag = 0; %default = unaligned data
end

%-----------------------------------------------------------------
%compute spatial filters based on covariance of the average
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
cfg.latency      = [-0.2 0.8-1./256];
%cfg.covariance   = 'yes';
%cfg.covariancewindow = [-0.2 0.8];
cfg.covariance   = 'no';
cfg.blc          = 'yes';
cfg.blcwindow    = [-0.2 0];
cfg.channel      = 'MEG';
tlck             = timelockanalysis(cfg, data);
clear data;

cfg.covariance       = 'yes';
cfg.covariancewindow = [-0.2 0.8-1./256];
cfg.vartrllength     = 'no';
tlck                 = timelockanalysis(cfg, tlck);

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
