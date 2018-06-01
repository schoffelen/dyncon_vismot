function [statl, statr, stat] = doSourceanalysisDICSpre2(subject, frequency, smoothing)

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

cd('/analyse/4/Project0030/roi');
%load([subject.name,'roimaskInducedGamma1']);
%load([subject.name,'roimaskTFCE']);
load([subject.name,'roimaskCongruency']);
sel1 = find(ismember(grid.inside,roi1));
sel2 = find(ismember(grid.inside,roi2));

indxl = [btrl(1):etrl(1) btrl(3):etrl(3)];
indxr = [btrl(4):etrl(4) btrl(2):etrl(2)];

%%median split
%rt    = freq.rt(indx);
%rt    = rt(:)';
%mrt   = median(rt);
%rt(rt<mrt)  = 1; %fast
%rt(rt>=mrt) = 2; %slow
%tmpsd = selectdata(sd, 'rpt', indx);

%get orientation in roi as consistent as possible to allow for complex averaging
sdorig = sd;
sd = source2sparse(sd);
ori1 = cat(2, sd.avg.ori{sel1});
ori2 = cat(2, sd.avg.ori{sel2});

sign1 = sum(repmat([1 2 4]', [1 size(ori1,2)]).*double(sign(ori1)==1));
sign2 = sum(repmat([1 2 4]', [1 size(ori2,2)]).*double(sign(ori2)==1));

m1 = mode(sign1);
m2 = mode(sign2);

flipori1 = zeros(1,size(ori1,2));
if m1==0,
  flipori1 = ismember(sign1, [3 5 6 7]);
elseif m1==1,
  flipori1 = ismember(sign1, [2 4 6 7]);
elseif m1==2,
  flipori1 = ismember(sign1, [1 4 5 7]);
elseif m1==3,
  flipori1 = ismember(sign1, [0 4 5 6]);
elseif m1==4,
  flipori1 = ismember(sign1, [1 2 3 7]);
elseif m1==5,
  flipori1 = ismember(sign1, [0 2 3 6]);
elseif m1==6,
  flipori1 = ismember(sign1, [0 1 3 5]);
elseif m1==7,
  flipori1 = ismember(sign1, [0 1 2 4]);
end
flipori2 = zeros(1,size(ori2,2));
if m2==0,
  flipori2 = ismember(sign2, [3 5 6 7]);
elseif m2==1,
  flipori2 = ismember(sign2, [2 4 6 7]);
elseif m2==2,
  flipori2 = ismember(sign2, [1 4 5 7]);
elseif m2==3,
  flipori2 = ismember(sign2, [0 4 5 6]);
elseif m2==4,
  flipori2 = ismember(sign2, [1 2 3 7]);
elseif m2==5,
  flipori2 = ismember(sign2, [0 2 3 6]);
elseif m2==6,
  flipori2 = ismember(sign2, [0 1 3 5]);
elseif m2==7,
  flipori2 = ismember(sign2, [0 1 2 4]);
end

flipori1 = -double(flipori1).*2+1;
flipori2 = -double(flipori2).*2+1;

for k = 1:3
  if k==1,
    indx = indxl;
  elseif k==2,
    indx = indxr;
  elseif k==3,
    indx = btrl(1):etrl(4);
  end

  if k<3,
    %exclude outliers >1.5 iqr from median
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
      indxl1 = indxl(sel(indx1));
      indxl2 = indxl(sel(indx2));
    elseif k==2,
      indxr1 = indxr(sel(indx1));
      indxr2 = indxr(sel(indx2));
    end
  else
    mrt1  = [];
    mrt2  = [];
    indx1 = [indxl1(:)' indxr1(:)'];
    indx2 = [indxl2(:)' indxr2(:)'];
    tmpsd = source;
  end
  tmpsd.cfg.trl = zeros(numel(tmpsd.cumtapcnt),3);
  tmpsd1 = selectdata(tmpsd, 'rpt', indx1);
  tmpsd2 = selectdata(tmpsd, 'rpt', indx2);
  
  cfgc = [];
  cfgc.method = 'coh';
  cfgc.refindx = sel1;
  
  scoh1 = ft_connectivityanalysis(cfgc, tmpsd1);
  scoh2 = ft_connectivityanalysis(cfgc, tmpsd2);
  cfgc.refindx = sel2;
  scoh3 = ft_connectivityanalysis(cfgc, tmpsd1);
  scoh4 = ft_connectivityanalysis(cfgc, tmpsd2);
  
  df1   = sum(tmpsd1.cumtapcnt(:))*2;
  df2   = sum(tmpsd2.cumtapcnt(:))*2;
  denom = sqrt(1./(df1-2) + 1./(df2-2));
  rsh   = [numel(scoh1.cohspctrm)./numel(sel1) numel(sel1)];
  %%dcoh1 = (atanh(trimmean(reshape(scoh1.cohspctrm, rsh),0.2,2)) - ...
  %%         atanh(trimmean(reshape(scoh2.cohspctrm, rsh),0.2,2)))./denom;
  dcoh1 = trimmean(abs(atanh(reshape(scoh1.cohspctrm, rsh))) - ...
           abs(atanh(reshape(scoh2.cohspctrm, rsh))),0.2,2)./denom;
  %dcoh1 = (atanh(abs(trimmean(repmat(flipori1, [rsh(1) 1]).*reshape(scoh1.cohspctrm, rsh),0.2,2))) - ...
  %         atanh(abs(trimmean(repmat(flipori1, [rsh(1) 1]).*reshape(scoh2.cohspctrm, rsh),0.2,2))))./denom;
  rsh   = [numel(scoh3.cohspctrm)./numel(sel2) numel(sel2)];
  %%dcoh2 = (atanh(trimmean(reshape(scoh3.cohspctrm, rsh),0.2,2)) - ...
  %%         atanh(trimmean(reshape(scoh4.cohspctrm, rsh),0.2,2)))./denom;
  dcoh2 = trimmean(abs(atanh(reshape(scoh3.cohspctrm, rsh))) - ...
          abs(atanh(reshape(scoh4.cohspctrm, rsh))),0.2,2)./denom;
  %dcoh2 = (atanh(abs(trimmean(repmat(flipori2, [rsh(1) 1]).*reshape(scoh3.cohspctrm, rsh),0.2,2))) - ...
  %         atanh(abs(trimmean(repmat(flipori2, [rsh(1) 1]).*reshape(scoh4.cohspctrm, rsh),0.2,2))))./denom;
  
  dcoh1(~isfinite(dcoh1)) = 0;
  dcoh2(~isfinite(dcoh2)) = 0; %account for nans

  dum = zeros(dim);
  dum(grid.inside) = dcoh1;
  dum(sdorig.outside)  = 0;
  dcoh1s           = smooth_vol(dum(:), krn, dim, sdorig.inside); 
  dcoh1            = dum(:);

  dum = zeros(dim);
  dum(grid.inside) = dcoh2;
  dum(sdorig.outside)  = 0;
  dcoh2s           = smooth_vol(dum(:), krn, dim, sdorig.inside); 
  dcoh2            = dum(:);

  stat         = [];
  stat.pos     = sdorig.pos;
  stat.inside  = sdorig.inside;
  stat.outside = sdorig.outside;
  stat.dim     = dim;
  stat.rt      = [mrt1 mrt2];
  stat.dcoh1   = dcoh1;
  stat.dcoh2   = dcoh2;
  stat.dcoh1s  = dcoh1s;
  stat.dcoh2s  = dcoh2s;
  stat.roi1    = roi1;
  stat.roi2    = roi2;
  stat.fwhm    = fwhm(:);
  
  if k==1,
    statl = stat;
  elseif k==2,
    statr = stat;
  elseif k==3,
    stat = stat;
  end
end
