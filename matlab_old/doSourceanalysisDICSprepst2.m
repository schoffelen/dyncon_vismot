function [stat, stat13, stat42, stat13b, stat42b, stat13x, stat42x] = doSourceanalysisDICSprepst2(subject, frequency, flag)

%Try to do the same as doSourceanalysisDICSprepst but then more memory efficient
%(accommodate for gamma-band range, high number of tapers)
%
%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's baseline and post
%stimulus interval are matched, and number of samples across conditions which are contrasted are matched)
%spatial filters are computed on all trials concatenated
%output:
%stat: act vs baseline depsamplesT
%stat13 stat42: act congruency effect
%stat13b stat42b: baselines congruency effect
%stat13x stat42x: relative change (log10(act)-log10(bas)) congruency effect

if nargin<3,
  flag = 0;
end

cd([subject.pathname,'freq/']);
if flag==0 && frequency<37,
  load([subject.name,'mtmfft004b']);
  %load([subject.name,'hanning'], 'freqpst');
elseif flag==0,
  load([subject.name,'mtmfft012b']);
elseif flag==1 && frequency<37,
  load([subject.name,'mtmfft_aligned004b']);
elseif flag==2 && frequency<37,
  load([subject.name,'mtmfft_aligned004bc']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned012b']);
elseif flag==2,
  load([subject.name,'mtmfft_aligned012bc']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
for k = 1:5
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+5) = length(freqpst{k}.cumtapcnt);
end
oklabel = gradlabel(ok==5);
for k = 1:5
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', frequency);
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

freq{1} = selectdata(freqpre{1:5}, 'param', 'fourierspctrm');
freq{2} = selectdata(freqpst{1:5}, 'param', 'fourierspctrm');
freq    = selectdata(freq{:}, 'param', 'fourierspctrm');
clear freqpre freqpst

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
cfg.lambda     = '5%';
cfg.feedback   = 'none';
%[u,s,v] = svd(tlckcovpre.cov);
%cfg.subspace   = pinv(sqrtm(tlckcovpre.cov+eye(length(tlckcovpre.label))*0.01*s(1,1)));
%cfg.lambda     = 1;
source         = sourceanalysis(cfg, freq);

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);
%1:5 are the baselines for conditions 1:5
%6:10 are the post stim intervals for conditions 1:5

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
krn  = compute_kernel(sd, 'truncate', 5e-5);
sd   = rmfield(sd, 'fwhm');
sd   = rmfield(sd, 'leadfield');
sd   = rmfield(sd, 'rough');
sd   = rmfield(sd, 'q');

%keep track of this
outside = sd.outside;
inside  = sd.inside;
nchunks = 6;
div     = floor(linspace(0,length(inside),nchunks+1));

cfg.method      = 'pcc';
cfg.keepmom     = 'yes';

%then divide the dipole grid into chunks, to be glued together afterwards
for kk = 1:nchunks  
  cfg.grid.inside = inside((div(kk)+1):div(kk+1));
  %project the single trials through the filters
  tmpsource       = sourceanalysis(cfg, freq);
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

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.statistic = 'indepsamplesT';
cfg2.parameter = 'pow';

%activation versus baseline on all trials
cfg2.statistic = 'depsamplesT';
cfg2.design    = ones(1,size(sd.pow,1));
cfg2.design(1:etrl(5)) = 2;
cfg2.design(2,:) = [1:etrl(5) 1:etrl(5)];
cfg2.ivar      = 1;
cfg2.uvar      = 2;
stat          = sourcestatistics(cfg2, sd);
stat.inside   = inside;
stat.outside  = outside;
stat.pos      = grid.pos;
dum           = zeros(size(stat.pos,1),1);
dum(inside)   = stat.stat;
stat.stat     = dum;
stat          = rmfield(stat,'mask');
stat          = rmfield(stat,'prob');
stat          = rmfield(stat,'dim');
dim           = pos2dim3d(stat.pos);
stat.stat2    = smooth_vol(stat.stat,krn,dim,inside);

%activation condition 1 vs 3
cfg2.statistic = 'indepsamplesT';
cfg2           = rmfield(cfg2,'uvar');
tmpsd = selectdata(sd, 'rpt', [btrl(6):etrl(6) btrl(8):etrl(8)]);
n1    = numel(btrl(6):etrl(6));
n2    = numel(btrl(8):etrl(8));
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat13     = sourcestatistics(cfg2, tmpsd);
stat13.inside  = inside;
stat13.outside = outside;
stat13.pos     = grid.pos;
dum           = zeros(size(stat13.pos,1),1);
dum(inside)   = stat13.stat;
stat13.stat   = dum;
stat13        = rmfield(stat13,'mask');
stat13        = rmfield(stat13,'prob');
stat13        = rmfield(stat13,'dim');
stat13.stat2  = smooth_vol(stat13.stat,krn,dim,inside);
stat13.cfg2   = [];

%activation condition 4 vs 2
tmpsd = selectdata(sd, 'rpt', [btrl(9):etrl(9) btrl(7):etrl(7)]);
n1    = numel(btrl(9):etrl(9));
n2    = numel(btrl(7):etrl(7));
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat42     = sourcestatistics(cfg2, tmpsd);
stat42.inside  = inside;
stat42.outside = outside;
stat42.pos     = grid.pos;
dum           = zeros(size(stat42.pos,1),1);
dum(inside)   = stat42.stat;
stat42.stat   = dum;
stat42        = rmfield(stat42,'mask');
stat42        = rmfield(stat42,'prob');
stat42        = rmfield(stat42,'dim');
stat42.stat2 = smooth_vol(stat42.stat,krn,dim,inside);
stat42.cfg2   = [];

%bsl condition 1 vs 3
tmpsd = selectdata(sd, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
n1    = numel(btrl(1):etrl(1));
n2    = numel(btrl(3):etrl(3));
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat13b    = sourcestatistics(cfg2, tmpsd);
stat13b.inside  = inside;
stat13b.outside = outside;
stat13b.pos     = grid.pos;
dum           = zeros(size(stat13b.pos,1),1);
dum(inside)   = stat13b.stat;
stat13b.stat   = dum;
stat13b        = rmfield(stat13b,'mask');
stat13b        = rmfield(stat13b,'prob');
stat13b        = rmfield(stat13b,'dim');
stat13b.stat2 = smooth_vol(stat13b.stat,krn,dim,inside);
stat13b.cfg2   = [];

%bsl condition 4 vs 2
tmpsd = selectdata(sd, 'rpt', [btrl(4):etrl(4) btrl(2):etrl(2)]);
n1    = numel(btrl(4):etrl(4));
n2    = numel(btrl(2):etrl(2));
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat42b    = sourcestatistics(cfg2, tmpsd);
stat42b.inside  = inside;
stat42b.outside = outside;
stat42b.pos     = grid.pos;
dum           = zeros(size(stat42b.pos,1),1);
dum(inside)   = stat42b.stat;
stat42b.stat   = dum;
stat42b        = rmfield(stat42b,'mask');
stat42b        = rmfield(stat42b,'prob');
stat42b        = rmfield(stat42b,'dim');
stat42b.stat2 = smooth_vol(stat42b.stat,krn,dim,inside);
stat42b.cfg2   = [];

%relative change condition 1 vs 3
tmpsd  = selectdata(sd, 'rpt', [btrl(6):etrl(6) btrl(8):etrl(8)]);
tmpsdb = selectdata(sd, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
n1     = numel(btrl(6):etrl(6));
n2     = numel(btrl(8):etrl(8));
tmpsd.pow  = log10(tmpsd.pow)-log10(tmpsdb.pow);
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat13x = sourcestatistics(cfg2, tmpsd);
stat13x.inside  = inside;
stat13x.outside = outside;
stat13x.pos     = grid.pos;
dum           = zeros(size(stat13x.pos,1),1);
dum(inside)   = stat13x.stat;
stat13x.stat   = dum;
stat13x        = rmfield(stat13x,'mask');
stat13x        = rmfield(stat13x,'prob');
stat13x        = rmfield(stat13x,'dim');
stat13x.stat2 = smooth_vol(stat13x.stat,krn,dim,inside);
stat13x.cfg2   = [];

%relative change condition 4 vs 2
tmpsd  = selectdata(sd, 'rpt', [btrl(9):etrl(9) btrl(7):etrl(7)]);
tmpsdb = selectdata(sd, 'rpt', [btrl(4):etrl(4) btrl(2):etrl(2)]);
n1     = numel(btrl(9):etrl(9));
n2     = numel(btrl(7):etrl(7));
tmpsd.pow  = log10(tmpsd.pow)-log10(tmpsdb.pow);
cfg2.design = [ones(1,n1) ones(1,n2)*2];
stat42x = sourcestatistics(cfg2, tmpsd);
stat42x.inside  = inside;
stat42x.outside = outside;
stat42x.pos     = grid.pos;
dum           = zeros(size(stat42x.pos,1),1);
dum(inside)   = stat42x.stat;
stat42x.stat   = dum;
stat42x        = rmfield(stat42x,'mask');
stat42x        = rmfield(stat42x,'prob');
stat42x        = rmfield(stat42x,'dim');
stat42x.stat2 = smooth_vol(stat42x.stat,krn,dim,inside);
stat42x.cfg2   = [];
