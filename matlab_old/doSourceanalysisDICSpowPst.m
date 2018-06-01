function [statl, statr] = doSourceanalysisDICSpowPst(subject, frequency)

%do sourceanalysis using all pst stimulus data 

%output:
%statl statr, congruency effect response left, congruency effect response right


fieldtripdefs

addpath /home/jan/matlab/robuststats

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/';

cd([datadir,'/freq/']);
if 1,
  if frequency>=40
    load([subject.name,'freqPst500']);
  else
    %FIXME compute 500 ms for low frequencies as well
    load([subject.name,'freqPst500low']);
  end


  %concatenate data and keep track of original trial numbers
  %gradlabel = freqpst{1}.grad.label;
  %ok = zeros(248,1);
  rtall = [];
  for k = 1:4
  %  [a,b] = match_str(gradlabel, freqpst{k}.label);
  %  ok(a) = ok(a)+1;
    ntrl(1,k) = length(freqpst{k}.cumtapcnt);
  end
  nmin = min(ntrl);
  for k = 1:4
    sel   = randperm(ntrl(k));
    sel   = sel(1:nmin);
    freqpst{k} = ft_selectdata(freqpst{k}, 'rpt', sort(sel));
    rt    = freqpst{k}.trialinfo(:,3);
    rtall = [rtall; rt];
    freqpst{k}.rt = rt./fsample;
  end
  
  %oklabel = gradlabel(ok==4);
  for k = 1:4
    warning off
    freqpst{k} = ft_selectdata(struct2double(freqpst{k}), 'foilim', frequency);
    warning on;
  %  freqpst{k} = ft_selectdata(freqpst{k}, 'channel', oklabel);
  end
  
  freq    = ft_selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
  freq.rt = rtall;
  clear freqpst
  ntrl(:) = nmin;
else
  load([subject.name,'freqPst400clean']);
  for k = 1:4
    ntrl(1,k) = sum(freq.trialinfo(:,1)==k);
  end
  rtall = freq.trialinfo(:,3);
  rt    = rtall;
  freq.rt = rt./fsample;
end

load([subject.pathname,'data/',subject.name,'data.mat'],'data1');
grad = data1.grad;
clear data1;

%load grid and prune leadfields
[a,b] = match_str(freq.label, grad.label);
load([subject.pathname,'grid/',subject.name,'grid6mm.mat']);
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
%freq.grad = convert_units(freq.grad, 'cm');

freq.grad = convert_units(grad, 'cm');

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
source         = ft_sourceanalysis(cfg, freq);

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
%krn  = compute_kernel(sd, 'truncate', 1e-5);
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
source          = ft_sourceanalysis(cfg, freq);
source          = checkdata(source, 'sourcerepresentation', 'new');
%source          = source2sparse(source);

%cd('/analyse/4/Project0030/roi');
%load([subject.name,'roimaskInducedGamma1']);
%load([subject.name,'roimaskTFCE']);
%load([subject.name,'roimaskCongruency']);
%sel1 = find(ismember(grid.inside,roi1));
%sel2 = find(ismember(grid.inside,roi2));

indxlc = [btrl(1):etrl(1)];
indxli = [btrl(3):etrl(3)];
indxrc = [btrl(4):etrl(4)];
indxri = [btrl(2):etrl(2)];

%get orientation in roi as consistent as possible to allow for complex averaging
sdorig = sd;
%sd = source2sparse(sd);
%ori1 = cat(2, sd.avg.ori{sel1});
%ori2 = cat(2, sd.avg.ori{sel2});

for k = 1:2
  if k==1,
    indx1 = indxlc;
    indx2 = indxli;
  elseif k==2,
    indx1 = indxrc;
    indx2 = indxri;
  end

  sd.cfg.trl = zeros(numel(sd.cumtapcnt),3);
  tmpsd1 = ft_checkdata(ft_selectdata(source, 'rpt', indx1), 'sourcerepresentation', 'new', 'haspow', 'yes');
  tmpsd2 = ft_checkdata(ft_selectdata(source, 'rpt', indx2), 'sourcerepresentation', 'new', 'haspow', 'yes');
  dat{1} = tmpsd1; n1 = numel(dat{1}.cumtapcnt);
  dat{2} = tmpsd2; n2 = numel(dat{2}.cumtapcnt);
  dat    = ft_selectdata(dat{:}, 'param', 'pow');
  df1    = sum(tmpsd1.cumtapcnt(:));
  df2    = sum(tmpsd2.cumtapcnt(:));
  clear tmpsd1 tmpsd2;
  
  cfgc                  = [];
  cfgc.implementation   = 'new';
  cfgc.statistic        = 'indepsamplesT';
  cfgc.method           = 'montecarlo';
  cfgc.parameter        = 'pow';
  cfgc.numrandomization = 0;
  cfgc.design           = [ones(1,n1) ones(1,n2)*2];
  stat                  = ft_sourcestatistics(cfgc, dat);
  clear dat;

  stat.inside  = sdorig.inside;
  stat.outside = sdorig.outside;
  stat.dim     = dim;
  stat.fwhm    = fwhm(:);
  stat.df1     = df1;
  stat.df2     = df2;

  if k==1,
    statl = stat;
  elseif k==2,
    statr = stat;
  end
end
