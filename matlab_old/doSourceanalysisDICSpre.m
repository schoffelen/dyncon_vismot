function [stat13, stat42] = doSourceanalysisDICSpre(subject, frequency, smoothing)

%do sourceanalysis using all pre stimulus data 

%output:
%Median split analysis on reaction time
%stat13 stat42;

fieldtripdefs

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

cd([subject.pathname,'freq/']);
if smoothing==0
  load([subject.name,'mtmfftPre000']);
elseif smoothing==3.75
  load([subject.name,'mtmfftPre038']);
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
  rt    = trl(:,4);
  rtall = [rtall; rt];
  freqpre{k}.rt = rt./fsample;
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
freq.rt = rtall;
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
krn  = compute_kernel(sd, 'truncate', 5e-5);
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
cfg2.statistic = 'indepsamplesT';
cfg2.ivar      = 1;

cfg3           = cfg2;
cfg3.statistic = 'yuent';

%congruency T-test
%activation condition 1 vs 3 response left
%activation condition 4 vs 2 response right
conditions = [1 3;4 2];
for mm = 1:2
    
  indx1 = btrl(conditions(mm,1)):etrl(conditions(mm,1));
  indx2 = btrl(conditions(mm,2)):etrl(conditions(mm,2));
  indx  = [indx1 indx2];
  
  %%median split
  %rt    = freq.rt(indx);
  %rt    = rt(:)';
  %mrt   = median(rt);
  %rt(rt<mrt)  = 1; %fast
  %rt(rt>=mrt) = 2; %slow
  %tmpsd = selectdata(sd, 'rpt', indx);
  
  %exclude outliers >1.5 iqr from median
  tmpsd = selectdata(sd, 'rpt', indx);
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
  nx    = ceil(numel(tmprt)/3);
  indx  = srtix([1:nx (numel(srtix)-nx+1):numel(srtix)]);
  tmprt = tmprt(indx);
  tmpsd = selectdata(tmpsd, 'rpt', indx);
  rt    = [ones(1,nx) ones(1,nx)*2];

  tmpsd.pow = standardise(log10(tmpsd.pow), 1);
  cfg2.design      = rt; 
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
  %tmpstat         = rmfield(tmpstat,'dim');
  tmpstat.stat2   = smooth_vol(tmpstat.stat,krn,dim,inside);
  tmpstat.cfg     = [];
  tmpstaty        = sourcestatistics(cfg3, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.staty   = dum;
  tmpstat.stat2y  = smooth_vol(tmpstat.staty,krn,dim,inside); 
  
  if mm==1,
    stat13 = tmpstat;
  elseif mm==2,
    stat42 = tmpstat;
  end
end
