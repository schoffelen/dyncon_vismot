function [stat13, stat42, stat13b, stat42b] = doSourceanalysisDICSprepst4(subject, frequency, flag)

%Try to do the something similar as doSourceanalysisDICSprepst2
%
%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's baseline and post
%stimulus interval are matched, and number of samples across conditions which are contrasted are matched)
%spatial filters are computed on all trials concatenated
%output:
%GLM analysis with only reaction time put in two regressors
%after mean correction of dependent variables 
%essentially this tests for an interaction between condition
%and rt, or in other words it tests the difference in correlation
%between power and reaction time across conditions
%stat13 stat42: act congruency effect in field stat and stat2
%stat13 stat42: act rt effect in field statrt stat2rt
%stat13b stat42b: as in act but then baselines

if nargin<3,
  flag = 0;
end

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

if length(subject.runnames)>1,
  for k = 1:length(subject.runnames)
    fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{k},subject.datafile];
    event{k}.event = read_event(fname);
  end
else
  event = read_event(fname);
end

cd([subject.pathname,'freq/']);
if flag==0 && frequency<37,
  load([subject.name,'mtmfft004bnew']);
elseif flag==0,
  load([subject.name,'mtmfft004bnew']);
  for k = 1:5
    trl1{k,1} = findcfg(freqpre{k}.cfg, 'trl');
    trl2{k,1} = findcfg(freqpst{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl1{k};
    freqpst{k}.cfg.trl = trl2{k};
  end
elseif flag==1 && frequency<37,
  load([subject.name,'mtmfft_aligned004bnew']);
elseif flag==2 && frequency<37,
  load([subject.name,'mtmfft_aligned004bcnew']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned004bnew']);
  for k = 1:5
    trl1{k,1} = findcfg(freqpre{k}.cfg, 'trl');
    trl2{k,1} = findcfg(freqpst{k}.cfg, 'trl');
  end
  load([subject.name,'mtmfft_aligned012b']);
  for k = 1:5
    freqpre{k}.cfg.trl = trl1{k};
    freqpst{k}.cfg.trl = trl2{k};
  end
elseif flag==2,
  load([subject.name,'mtmfft_aligned012bcnew']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
rtall = [];
for k = 1:5
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+5) = length(freqpst{k}.cumtapcnt);
  if k==5.
    freqpre{k}.rt = nan+zeros(length(freqpre{k}.cumtapcnt),1);
    freqpst{k}.rt = nan+zeros(length(freqpst{k}.cumtapcnt),1);
  else
    rt = data2rt(freqpre{k}, event, 2);
    freqpre{k}.rt = rt./fsample;
    freqpst{k}.rt = rt./fsample;
  end
  rtall = [rtall;rt./fsample]; %rt for condition 5 does not mean anything
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
freq.rt = [rtall;rtall];
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
dim     = pos2dim3d(newgrid.pos);
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
cfg2.parameter = 'pow';
cfg2.statistic = 'glm';
cfg2.glm.statistic = 'T';

%congruency and RT glm for response left
%activation condition 1 vs 3
n1        = numel(btrl(6):etrl(6));
n2        = numel(btrl(8):etrl(8));
tmpsd     = selectdata(sd, 'rpt', [btrl(6):etrl(6) btrl(8):etrl(8)]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
m1 = mean(tmpsd.pow(1:n1,:));
m2 = mean(tmpsd.pow((n1+1):(n1+n2),:));
tmpsd.pow(1:n1,:)           = tmpsd.pow(1:n1,:)           - m1(ones(1,n1),:);
tmpsd.pow((n1+1):(n1+n2),:) = tmpsd.pow((n1+1):(n1+n2),:) - m2(ones(1,n2),:);
cfg2.design                   = zeros(2,n1+n2);
cfg2.design(1,1:n1)           = freq.rt(btrl(6):etrl(6))';
cfg2.design(2,(n1+1):(n1+n2)) = freq.rt(btrl(8):etrl(8))';
cfg2.glm.contrast             = [1 -1];
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
stat13.cfg   = [];

%congruency and RT glm for response right
%activation condition 4 vs 2
n1        = numel(btrl(9):etrl(9));
n2        = numel(btrl(7):etrl(7));
tmpsd     = selectdata(sd, 'rpt', [btrl(9):etrl(9) btrl(7):etrl(7)]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
m1 = mean(tmpsd.pow(1:n1,:));
m2 = mean(tmpsd.pow((n1+1):(n1+n2),:));
tmpsd.pow(1:n1,:)           = tmpsd.pow(1:n1,:) - m1(ones(1,n1),:);
tmpsd.pow((n1+1):(n1+n2),:) = tmpsd.pow((n1+1):(n1+n2),:) - m2(ones(1,n2),:);
cfg2.design                   = zeros(2,n1+n2);
cfg2.design(1,1:n1)           = freq.rt(btrl(9):etrl(9))';
cfg2.design(2,(n1+1):(n1+n2)) = freq.rt(btrl(7):etrl(7))';
cfg2.glm.contrast             = [1 -1];
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
stat42.stat2  = smooth_vol(stat42.stat,krn,dim,inside);
stat42.cfg   = [];

%congruency and RT glm for response left baseline
%baseline condition 1 vs 3
n1        = numel(btrl(1):etrl(1));
n2        = numel(btrl(3):etrl(3));
tmpsd     = selectdata(sd, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
m1 = mean(tmpsd.pow(1:n1,:));
m2 = mean(tmpsd.pow((n1+1):(n1+n2),:));
tmpsd.pow(1:n1,:)           = tmpsd.pow(1:n1,:) - m1(ones(1,n1),:);
tmpsd.pow((n1+1):(n1+n2),:) = tmpsd.pow((n1+1):(n1+n2),:) - m2(ones(1,n2),:);
cfg2.design                   = zeros(2,n1+n2);
cfg2.design(1,1:n1)           = freq.rt(btrl(1):etrl(1))';
cfg2.design(2,(n1+1):(n1+n2)) = freq.rt(btrl(3):etrl(3))';
cfg2.glm.contrast  = [1 -1];
stat13b     = sourcestatistics(cfg2, tmpsd);
stat13b.inside  = inside;
stat13b.outside = outside;
stat13b.pos     = grid.pos;
dum           = zeros(size(stat13b.pos,1),1);
dum(inside)   = stat13b.stat;
stat13b.stat   = dum;
stat13b        = rmfield(stat13b,'mask');
stat13b        = rmfield(stat13b,'prob');
stat13b        = rmfield(stat13b,'dim');
stat13b.stat2  = smooth_vol(stat13b.stat,krn,dim,inside);
stat13b.cfg   = [];

%congruency and RT glm for response right baseline
%baseline condition 4 vs 2
n1        = numel(btrl(4):etrl(4));
n2        = numel(btrl(2):etrl(2));
tmpsd     = selectdata(sd, 'rpt', [btrl(4):etrl(4) btrl(2):etrl(2)]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
m1 = mean(tmpsd.pow(1:n1,:));
tmpsd.pow(1:n1,:) = tmpsd.pow(1:n1,:) - m1(ones(1,n1),:);
m2 = mean(tmpsd.pow((n1+1):(n1+n2),:));
tmpsd.pow((n1+1):(n1+n2),:) = tmpsd.pow((n1+1):(n1+n2),:) - m2(ones(1,n2),:);
cfg2.design                   = zeros(2,n1+n2);
cfg2.design(1,1:n1)           = freq.rt(btrl(4):etrl(4))';
cfg2.design(2,(n1+1):(n1+n2)) = freq.rt(btrl(2):etrl(2))';
cfg2.glm.contrast             = [1 -1];
stat42b     = sourcestatistics(cfg2, tmpsd);
stat42b.inside  = inside;
stat42b.outside = outside;
stat42b.pos     = grid.pos;
dum           = zeros(size(stat42b.pos,1),1);
dum(inside)   = stat42b.stat;
stat42b.stat   = dum;
stat42b        = rmfield(stat42b,'mask');
stat42b        = rmfield(stat42b,'prob');
stat42b        = rmfield(stat42b,'dim');
stat42b.stat2  = smooth_vol(stat42b.stat,krn,dim,inside);
stat42b.cfg    = [];
