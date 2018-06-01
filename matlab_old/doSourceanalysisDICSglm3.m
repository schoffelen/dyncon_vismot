function [stat13, stat42, stat13b, stat42b] = doSourceanalysisDICSglm3(subject, frequency)

%Try to do the something similar as doSourceanalysisDICSglm2
%
%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's %baseline and post stimulus interval are matched, and number of samples across conditions which are
%contrasted are matched) spatial filters are computed on all trials concatenated. overall trials are
%stratified for reaction time
%output:
%GLM analysis
%stat13 stat42: act congruency effect in field stat and stat2
%stat13 stat42: act rt effect in field statrt stat2rt
%stat13b stat42b: as in act but then baselines

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
blockpre = [];
blockpst = [];
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+4) = length(freqpst{k}.cumtapcnt);
  trl   = findcfg(freqpre{k}.cfg, 'trl');
  rt    = trl(:,4);
  freqpre{k}.rt = rt./fsample;
  blockpre = [blockpre trl2blockregressor(subject, trl)];
  
  trl   = findcfg(freqpst{k}.cfg, 'trl');
  rt    = trl(:,4);
  rtall = [rtall; rt];
  freqpst{k}.rt = rt./fsample;
  blockpst = [blockpst trl2blockregressor(subject, trl)];
end
block = [blockpre blockpst];

oklabel = gradlabel(ok==4);
for k = 1:4
  warning off
  freqpst{k} = selectdata(struct2double(freqpst{k}), 'foilim', frequency);
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

freq{1} = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
freq{2} = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
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

%congruency and RT glm
%activation condition 1 vs 3 response left
%activation condition 4 vs 2 response right
%baseline condition 1 vs 3 response left
%baseline condition 4 vs 2 response right
conditions = [5 7;8 6;1 3;4 2];
for mm = 1:4
    
  indx1 = btrl(conditions(mm,1)):etrl(conditions(mm,1));
  indx2 = btrl(conditions(mm,2)):etrl(conditions(mm,2));
  indx  = [indx1 indx2];
  tmpsd = selectdata(sd, 'rpt', indx);
  tmpsd.pow = standardise(log10(tmpsd.pow), 1);
  n1    = numel(indx1);
  n2    = numel(indx2);
  cfg2.design = block(:, indx);
  cfg2.design(end+1,:) = freq.cumtapcnt(indx)';
  %cfg2.design(end+1,:) = (freq.cumtapcnt(indx)').^2;
  cfg2.design(end+1,:) = [ones(1,n1) -ones(1,n2)]; %congruency regressor
  cfg2.design(end+1,:) = freq.rt(indx)';   %rt regressor
  cfg2.design          = standardise(orthogonalise(cfg2.design')',2);
  cfg2.glm.contrast  = zeros(1,size(cfg2.design,1));
  cfg2.glm.contrast(end-1:end) = [1 0]; %congruency contrast
  
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
  cfg2.glm.contrast(end-1:end) = [0 1]; %correlation with rt
  tmpstatrt       = sourcestatistics(cfg2, tmpsd);
  dum(inside)     = tmpstatrt.stat;
  tmpstat.statrt  = dum;
  tmpstat.stat2rt = smooth_vol(tmpstat.statrt,krn,dim,inside); 
  
  if mm==1,
    stat13 = tmpstat;
  elseif mm==2,
    stat42 = tmpstat;
  elseif mm==3,
    stat13b = tmpstat;
  elseif mm==4,
    stat42b = tmpstat;
  end
end

function [x] = orthogonalise(input)

input = standardise(input,1);
x     = input(:,1);

for i = 2:size(input,2)
  y   = input(:,i);
  y   = y-x*(pinv(x'*x)*x'*y);
  if any(y), x = [x y]; end
end
