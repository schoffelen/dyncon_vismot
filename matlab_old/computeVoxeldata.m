function [sd, btrl, etrl, rtall] = computeVoxeldata(subject, frequency, smoothing)

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
elseif smoothing==5
  load([subject.name,'mtmfftPre400_050']);
elseif smoothing==3.2
  cd('/analyse/1/Project0002/tmpProject0030/freq/');
  load([subject.name,'mtmfftPre625_032']);
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

freq.grad = convert_units(freq.grad, 'cm');

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
  %newgrid.leadfield{indx} = tmplf(b,:);
  newgrid.leadfield{indx} = tmplf(b,1); %will not be used anyway; is necessary to survive
  %beamformer_pcc
end

%load([subject.pathname,'filter/',subject.name,'filt',num2str(round(10*smoothing),'%03d'), ...
%  '_',num2str(frequency,'%03d')]);
%load([subject.pathname,'filter/',subject.name,'filt400_',num2str(round(10*smoothing),'%03d'), ...
%  '_',num2str(frequency,'%03d')]);
load(['/analyse/1/Project0002/tmpProject0030/filter/',subject.name,'filt625_',num2str(round(10*smoothing),'%03d'), ...
  '_',num2str(round(10*frequency),'%03d')]);
eval('newgrid.filter = filt;');

vol       = read_vol([subject.pathname,'vol/',subject.name,'vol.mat']);
if strcmp(subject.datafile(1),'h'),
  %get vol in dewar space
  load([subject.pathname,'dewar2head_avg/',subject.name,'dewar2head_avg.mat'])
  vol       = transform_vol(inv(M), vol);
end

cfg            = [];
cfg.method     = 'pcc';
cfg.frequency  = frequency;
cfg.grid       = newgrid;
cfg.vol        = vol;
cfg.keepmom    = 'yes';
cfg.powmethod  = 'lambda1';
cfg.lambda     = '10%';
cfg.feedback   = 'none';

%adjust inside and outside
outside = find(~isfinite(fwhm));
inside  = find( isfinite(fwhm));
dim     = pos2dim3d(newgrid.pos);
nchunks = 4;
div     = floor(linspace(0,length(inside),nchunks+1));

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

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);
