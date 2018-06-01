function [stat13, stat42, stat13b, stat42b, fwhm] = doSourceanalysisDICSprepstStratified(subject, frequency, flag)

%Try to do the something similar as doSourceanalysisDICSprepst3
%
%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's baseline and post
%stimulus interval are matched, and number of samples across conditions which are contrasted are matched)
%spatial filters are computed on all trials concatenated
%output:
%GLM analysis
%based on data stratified for RT

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
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+4) = length(freqpst{k}.cumtapcnt);
  if k==5.
    rt = nan+zeros(length(freqpre{k}.cumtapcnt),1);
    freqpre{k}.rt = rt;
    rt = nan+zeros(length(freqpst{k}.cumtapcnt),1);
    freqpst{k}.rt = rt;
  else
    rt = data2rt(freqpre{k}, event, 2);
    freqpre{k}.rt = rt./fsample;
    freqpst{k}.rt = rt./fsample;
  end
  rtall = [rtall;rt./fsample]; %rt for condition 5 does not mean anything
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

load([subject.pathname,'stratifyRT/',subject.name,'stratifyRT.mat']);
for k = 1:4
  freqpre{k} = selectdata(freqpre{k}, 'rpt', find(isfinite(output{k})));
  freqpst{k} = selectdata(freqpst{k}, 'rpt', find(isfinite(output{k})));
  freqpre{k}.rt = freqpre{k}.rt(find(isfinite(output{k})));
  freqpst{k}.rt = freqpst{k}.rt(find(isfinite(output{k})));
  ntrl(k) = sum(isfinite(output{k}));
  ntrl(k+4) = sum(isfinite(output{k}));
end
freq{1} = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
freq{2} = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
freq    = selectdata(freq{:}, 'param', 'fourierspctrm');
rtall   = cat(2,output{:});
rtall   = rtall(isfinite(rtall));
freq.rt = log([rtall(:);rtall(:)]);
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
%1:4 are the baselines for conditions 1:4
%5:8 are the post stim intervals for conditions 1:4

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

%stratify for RT
sel1   = btrl(5):etrl(5);
sel2   = btrl(7):etrl(7);

sd.cfg.trl = zeros(length(freq.cumtapcnt),3);
tmpsd = selectdata(sd, 'rpt', [sel1 sel2]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
n1    = numel(sel1);
n2    = numel(sel2);
%orthogonalise with respect to the condition
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
rt1   = rt1(:)' - mean(rt1(:));
rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(2,:) = [rt1 rt2];
cfg2.glm.contrast  = [1 0];
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
cfg2.glm.contrast = [0 1]; %correlation with rt
stat13rt       = sourcestatistics(cfg2, tmpsd);
dum(inside)    = stat13rt.stat;
stat13.statrt  = dum;
stat13.stat2rt = smooth_vol(stat13.statrt,krn,dim,inside); 

%congruency and RT glm for response right
%activation condition 4 vs 2

%stratify for RT
sel1   = btrl(8):etrl(8);
sel2   = btrl(6):etrl(6);

tmpsd = selectdata(sd, 'rpt', [sel1 sel2]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
n1    = numel(sel1);
n2    = numel(sel2);
%orthogonalise with respect to condition
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
rt1   = rt1(:)' - mean(rt1(:));
rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(2,:) = [rt1 rt2];
cfg2.glm.contrast  = [1 0];
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
cfg2.glm.contrast = [0 1]; %correlation with rt
stat42rt       = sourcestatistics(cfg2, tmpsd);
dum(inside)    = stat42rt.stat;
stat42.statrt  = dum;
stat42.stat2rt = smooth_vol(stat42.statrt,krn,dim,inside); 

%congruency and RT glm for response left baseline
%baseline condition 1 vs 3

%stratify for RT
sel1   = btrl(1):etrl(1);
sel2   = btrl(3):etrl(3);

tmpsd = selectdata(sd, 'rpt', [sel1 sel2]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
n1    = numel(sel1);
n2    = numel(sel2);
%orthogonalise with respect to the condition
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
rt1   = rt1(:)' - mean(rt1(:));
rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(2,:) = [rt1 rt2];
cfg2.glm.contrast  = [1 0];
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
cfg2.glm.contrast = [0 1]; %correlation with rt
stat13brt       = sourcestatistics(cfg2, tmpsd);
dum(inside)    = stat13brt.stat;
stat13b.statrt  = dum;
stat13b.stat2rt = smooth_vol(stat13b.statrt,krn,dim,inside); 

%congruency and RT glm for response right baseline
%baseline condition 4 vs 2

%stratify for RT
sel1   = btrl(4):etrl(4);
sel2   = btrl(2):etrl(2);

tmpsd = selectdata(sd, 'rpt', [sel1 sel2]);
tmpsd.pow = standardise(log10(tmpsd.pow), 1);
n1    = numel(sel1);
n2    = numel(sel2);
%orthogonalise with respect to the condition
rt1   = freq.rt(sel1);
rt2   = freq.rt(sel2);
rt1   = rt1(:)' - mean(rt1(:));
rt2   = rt2(:)' - mean(rt2(:));
cfg2.design      = [ones(1,n1) -ones(1,n2)]; %congruency regressor
cfg2.design(2,:) = [rt1 rt2];
cfg2.glm.contrast  = [1 0];
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
stat42b.cfg   = [];
cfg2.glm.contrast = [0 1]; %correlation with rt
stat42brt       = sourcestatistics(cfg2, tmpsd);
dum(inside)    = stat42brt.stat;
stat42b.statrt  = dum;
stat42b.stat2rt = smooth_vol(stat42b.statrt,krn,dim,inside); 
