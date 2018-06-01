function [statlc, statli, statrc, statri] = doSourceanalysisDICSacPre(subject, frequency)

%do sourceanalysis using all pre stimulus data 

%output:
%Median split analysis on reaction time
%stat13 stat42;

fieldtripdefs

addpath /home/jan/matlab/robuststats

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/';

cd([datadir,'/freq/']);
if 1,
  if frequency<40
    load([subject.name,'freqPre400']);
  else  
    load([subject.name,'freqPre400high']);
  end
  %concatenate data and keep track of original trial numbers
  %gradlabel = freqpre{1}.grad.label;
  %ok = zeros(248,1);
  rtall = [];
  for k = 1:4
  %  [a,b] = match_str(gradlabel, freqpre{k}.label);
  %  ok(a) = ok(a)+1;
    ntrl(1,k) = length(freqpre{k}.cumtapcnt);
    rt    = freqpre{k}.trialinfo(:,3);
    rtall = [rtall; rt];
    freqpre{k}.rt = rt./fsample;
  end

  %oklabel = gradlabel(ok==4);
  for k = 1:4
    warning off
    freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
    warning on;
  %  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
  end
  
  freq    = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
  freq.rt = rtall;
  clear freqpre
else
  load([subject.name,'freqPre400clean']);
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
source          = source2sparse(source);

cd('/analyse/4/Project0030/roi');
%load([subject.name,'roimaskInducedGamma1']);
%load([subject.name,'roimaskTFCE']);
load([subject.name,'roimaskCongruency']);
sel1 = find(ismember(grid.inside,roi1));
sel2 = find(ismember(grid.inside,roi2));

indxlc = [btrl(1):etrl(1)];
indxli = [btrl(3):etrl(3)];
indxrc = [btrl(4):etrl(4)];
indxri = [btrl(2):etrl(2)];

%get orientation in roi as consistent as possible to allow for complex averaging
sdorig = sd;
sd = source2sparse(sd);
ori1 = cat(2, sd.avg.ori{sel1});
ori2 = cat(2, sd.avg.ori{sel2});

for k = 1:4
  if k==1,
    indx = indxlc;
  elseif k==2,
    indx = indxli;
  elseif k==3,
    indx = indxrc;
  elseif k==4,
    indx = indxri;
  end

  if k<5,
    %exclude outliers >1.5 iqr from median
    source.cfg.trl = zeros(numel(source.cumtapcnt),3);
    tmpsd = selectdata(source, 'rpt', indx);
    tmprt = freq.rt(indx);
    tmprt = tmprt(:)';
    mrt   = median(tmprt);
    irt   = iqr(tmprt);
    thrlo = mrt - 1.5*irt;
    thrhi = mrt + 1.5*irt;
    tmprt(tmprt<thrlo)=nan;
    tmprt(tmprt>thrhi)=nan;
    sel   = find(isfinite(tmprt));
    tmprt = tmprt(sel);
    tmpsd = selectdata(tmpsd, 'rpt', sel);
    [srt,srtix] = sort(tmprt);
    %nx    = ceil(numel(tmprt)/3);
    nx    = ceil(numel(tmprt)/2);
    indx1  = srtix([1:nx]);
    mrt1   = mean(srt(1:nx));
    indx2  = srtix([(numel(srtix)-nx+1):numel(srtix)]);
    mrt2   = mean(srt(end-nx+[1:nx]));
    if k==1,
      %indx1 = indxlc(sel(indx1));
      %indx2 = indxlc(sel(indx2));
      selx  = sel2; %lc means response cue in right hemi
   elseif k==2,
      %indx1 = indxli(sel(indx1));
      %indx2 = indxli(sel(indx2));
      %selx  = sel2;
      selx  = sel1; %li means response cue in left hemi
    elseif k==3,
      %indx1 = indxrc(sel(indx1));
      %indx2 = indxrc(sel(indx2));
      selx  = sel1; %rc means response cue in left hemi
    elseif k==4,
      %indx1 = indxri(sel(indx1));
      %indx2 = indxri(sel(indx2));
      %selx  = sel1;
      selx  = sel2; %ri means response cue in right hemi
    end
  else
  end
  tmpsd.cfg.trl = zeros(numel(tmpsd.cumtapcnt),3);
  tmpsd1 = selectdata(tmpsd, 'rpt', indx1);
  tmpsd2 = selectdata(tmpsd, 'rpt', indx2);
  
  cfgc = [];
  cfgc.method = 'amplcorr';
  cfgc.refindx = sel1;
  scoh1 = ft_connectivityanalysis(cfgc, tmpsd1);
  scoh2 = ft_connectivityanalysis(cfgc, tmpsd2);
  cfgc.refindx = sel2;
  scoh3 = ft_connectivityanalysis(cfgc, tmpsd1);
  scoh4 = ft_connectivityanalysis(cfgc, tmpsd2);

  nref  = numel(sel1);
  tmpx1 = zeros(prod(dim), nref);
  tmpx2 = zeros(prod(dim), nref);
  tmpx1(grid.inside, :) = reshape(scoh1.amplcorrspctrm, numel(scoh1.amplcorrspctrm)/nref, nref);
  tmpx2(grid.inside, :) = reshape(scoh2.amplcorrspctrm, numel(scoh2.amplcorrspctrm)/nref, nref);
  nref  = numel(sel2);
  tmpx3 = zeros(prod(dim), nref);
  tmpx4 = zeros(prod(dim), nref);
  tmpx3(grid.inside, :) = reshape(scoh3.amplcorrspctrm, numel(scoh3.amplcorrspctrm)/nref, nref);
  tmpx4(grid.inside, :) = reshape(scoh4.amplcorrspctrm, numel(scoh4.amplcorrspctrm)/nref, nref);

  tmpx1(sdorig.outside,:) = nan;
  tmpx2(sdorig.outside,:) = nan;
  tmpx3(sdorig.outside,:) = nan;
  tmpx4(sdorig.outside,:) = nan;

  df1   = sum(tmpsd1.cumtapcnt(:))*2;
  df2   = sum(tmpsd2.cumtapcnt(:))*2;
  df3   = sum(tmpsd1.cumtapcnt(:))*2;
  df4   = sum(tmpsd2.cumtapcnt(:))*2;
  denom = sqrt(1./(df1-2) + 1./(df2-2));
  
  stat         = [];
  stat.pos     = sdorig.pos;
  stat.inside  = sdorig.inside;
  stat.outside = sdorig.outside;
  stat.dim     = dim;
  stat.rt      = [mrt1 mrt2];
  stat.roi1    = roi1;
  stat.roi2    = roi2;
  stat.fwhm    = fwhm(:);
  stat.coh1    = tmpx1;
  stat.coh2    = tmpx2;
  stat.coh3    = tmpx3;
  stat.coh4    = tmpx4;
  stat.dimord  = 'pos_pos';
  stat.df1     = df1;
  stat.df2     = df2;
  stat.df3     = df3;
  stat.df4     = df4;

  if k==1,
    statlc = stat;
  elseif k==2,
    statli = stat;
  elseif k==3,
    statrc = stat;
  elseif k==4,
    statri = stat;
  end
end
