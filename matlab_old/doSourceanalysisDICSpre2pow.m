function [statlc, statli, statrc, statri] = doSourceanalysisDICSpre2pow(subject, frequency, smoothing)

%do sourceanalysis using all pre stimulus data 

%output:
%Median split analysis on reaction time
%stat13 stat42;

fieldtripdefs

addpath /home/jan/matlab/robuststats

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

cd([subject.pathname,'freq/']);
if smoothing==0
  load([subject.name,'mtmfftPre000']);
elseif smoothing==3.75
  load([subject.name,'mtmfftPre038']);
elseif smoothing==5
  load([subject.name,'mtmfftPre400_050']);
elseif smoothing==7.5
  load([subject.name,'mtmfftPre075']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpre{1}.grad.label;
ok = zeros(248,1);
rtall = [];
for k = 1:4
  [a,b] = match_str(gradlabel, freqpre{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  trl   = findcfg(freqpre{k}.cfg, 'trl');
  if ~isempty(trl)
    rt    = trl(:,4);
  else
    rt = freqpre{k}.trialinfo;
  end
  rtall = [rtall; rt];
  freqpre{k}.rt = rt(:,1)./fsample;
  %blockpst = [blockpst trl2blockregressor(subject, trl)];
end

oklabel = gradlabel(ok==4);
for k = 1:4
  warning off
  freqpre{k} = selectdata(struct2double(freqpre{k}), 'foilim', frequency);
  warning on;
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

freq    = selectdata(freqpre{1:4}, 'param', 'fourierspctrm');
freq.rt = rtall(:,1);
clear freqpre

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
source          = ft_sourceanalysis(cfg, freq);
source          = checkdata(source, 'sourcerepresentation', 'new');
source          = source2sparse(source);

indxlc = [btrl(1):etrl(1)];
indxli = [btrl(3):etrl(3)];
indxrc = [btrl(4):etrl(4)];
indxri = [btrl(2):etrl(2)];

sdorig = sd;
sd     = source2sparse(sd);

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

  %exclude outliers >1.5 iqr from median
  source.trialinfo = freq.trialinfo;
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
  nx    = floor(numel(tmprt)/2);
  indx1  = srtix([1:nx]);
  mrt1   = mean(srt(1:nx));
  indx2  = srtix([(numel(srtix)-nx+1):numel(srtix)]);
  mrt2   = mean(srt(end-nx+[1:nx]));
  
  tmpsd.cfg.trl = zeros(numel(tmpsd.cumtapcnt),3);
  tmpsd1 = selectdata(tmpsd, 'rpt', indx1); tmpsd1 = checkdata(tmpsd1, 'sourcerepresentation', 'new', 'haspow', 'yes');
  tmpsd2 = selectdata(tmpsd, 'rpt', indx2); tmpsd2 = checkdata(tmpsd2, 'sourcerepresentation', 'new', 'haspow', 'yes');
  
  cfgc = [];
  cfgc.method = 'montecarlo';
  cfgc.numrandomization  = 100;
  cfgc.statistic = 'yuent';
  cfgc.parameter = 'pow';
  cfgc.implementation = 'new';
  cfgc.design = [ones(1,size(tmpsd1.pow,2)) 2*ones(1,size(tmpsd2.pow,2))];
  tmpsd1.pow = log10(tmpsd1.pow);
  tmpsd2.pow = log10(tmpsd2.pow);
  tmpsd1.pow = [tmpsd1.pow tmpsd2.pow];

  tmpstat = ft_sourcestatistics(cfgc, tmpsd1);

  dum = zeros(dim);
  dum(grid.inside)  = tmpstat.stat;
  dum(sdorig.outside) = 0;
  tmp = smooth_vol(dum(:), krn, dim, sdorig.inside); 
  
  stat         = [];
  stat.pos     = sdorig.pos;
  stat.inside  = sdorig.inside;
  stat.outside = sdorig.outside;
  stat.dim     = dim;
  stat.rt      = [mrt1 mrt2];
  stat.stat    = dum(:);
  stat.stat2   = tmp(:);
  stat.fwhm    = fwhm(:);
  dum(grid.inside) = tmpstat.prob;
  stat.prob    = dum(:);
  
  % add information with respect to the length of the baselines
  ix = ismember(subject.trlid, tmpsd1.trialinfo(:,2));
  stat.bsl1 = -subject.startfix(ix) + subject.trigsmp(ix);
  ix = ismember(subject.trlid, tmpsd2.trialinfo(:,2));
  stat.bsl2 = -subject.startfix(ix) + subject.trigsmp(ix);
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
