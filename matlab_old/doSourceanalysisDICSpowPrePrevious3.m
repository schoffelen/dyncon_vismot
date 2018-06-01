function [statl, statr] = doSourceanalysisDICSpowPrePrevious3(subject, frequency)

%do sourceanalysis using all pre stimulus data 

%output:
%statl, statr

ft_defaults

addpath /home/jan/matlab/robuststats

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/';

cd([datadir,'/freq/']);
if 1,
  if frequency>=40
    load([subject.name,'freqPre400high']);
    freq = ft_selectdata(freqpre{:}, 'param', 'fourierspctrm');
    clear freqpre;
  elseif frequency>=16
    %FIXME compute 500 ms for low frequencies as well
    load([subject.name,'freqPre400clean']);
  else
    load([subject.name,'freqPre400hanning']);
    freq = ft_selectdata(freqpre{:}, 'param', 'fourierspctrm');
    clear freqpre;
  end


  %concatenate data and keep track of original trial numbers
  %gradlabel = freqpre{1}.grad.label;
  %ok = zeros(248,1);
  rtall = [];
  for k = 1:4
  %  [a,b] = match_str(gradlabel, freqpre{k}.label);
  %  ok(a) = ok(a)+1;
    ntrl(1,k) = sum(freq.trialinfo(:,1)==k);
  end
  %nmin = min(ntrl);
  %for k = 1:4
  %  sel   = randperm(ntrl(k));
  %  sel   = sel(1:nmin);
  %  freqpre{k} = ft_selectdata(freqpre{k}, 'rpt', sort(sel));
  %  rt    = freqpre{k}.trialinfo(:,3);
  %  rtall = [rtall; rt];
  %  freqpre{k}.rt = rt./fsample;
  %end
  %ntrl(:) = nmin;
  
  %oklabel = gradlabel(ok==4);
  warning off
  freq = ft_selectdata(struct2double(freq), 'foilim', frequency);
  warning on;
  
  freq.rt = rtall;
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
cfg.lambda     = '5%';
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

indxlc = [btrl(1):etrl(1)];
indxli = [btrl(3):etrl(3)];
indxrc = [btrl(4):etrl(4)];
indxri = [btrl(2):etrl(2)];

%get orientation in roi as consistent as possible to allow for complex averaging
sdorig = sd;

for k = 1:2
  nxx = size(subject.trlid,1);
  if k==1,
    trlx = freq.trialinfo(indxlc,2);
    ix1  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[1 4]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix1  = subject.trlid([zeros(nxx,1)==1 ix1]);
    indx1a  = indxlc(ismember(trlx, ix1));
    ix2  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[2 3]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix2  = subject.trlid([zeros(nxx,1)==1 ix2]);
    indx2a  = indxlc(ismember(trlx, ix2));
    
    trlx = freq.trialinfo(indxli,2);
    ix1  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[1 4]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix1  = subject.trlid([zeros(nxx,1)==1 ix1]);
    indx1b  = indxli(ismember(trlx, ix1));
    ix2  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[2 3]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix2  = subject.trlid([zeros(nxx,1)==1 ix2]);
    indx2b  = indxli(ismember(trlx, ix2));
  
    indx1 = [indx1a(:);indx1b(:)];
    indx2 = [indx2a(:);indx2b(:)];

  elseif k==2,
    trlx = freq.trialinfo(indxrc,2);
    ix1  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[1 4]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix1  = subject.trlid([zeros(nxx,1)==1 ix1]);
    indx1a  = indxrc(ismember(trlx, ix1));
    ix2  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[2 3]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix2  = subject.trlid([zeros(nxx,1)==1 ix2]);
    indx2a  = indxrc(ismember(trlx, ix2));
    
    trlx = freq.trialinfo(indxri,2);
    ix1  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[1 4]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix1  = subject.trlid([zeros(nxx,1)==1 ix1]);
    indx1b  = indxri(ismember(trlx, ix1));
    ix2  = ismember(subject.trlid(:,2:end),trlx) & ismember(subject.trl(:,1:3),[2 3]) & ...
           subject.correct(:,2:end) & subject.correct(:,1:3);
    ix2  = subject.trlid([zeros(nxx,1)==1 ix2]);
    indx2b  = indxri(ismember(trlx, ix2));
  
    indx1 = [indx1a(:);indx1b(:)];
    indx2 = [indx2a(:);indx2b(:)];
  
  end

  sd.cfg.trl = zeros(numel(sd.cumtapcnt),3);
  tmpsd1 = ft_checkdata(ft_selectdata(source, 'rpt', indx1a), 'sourcerepresentation', 'new', 'haspow', 'yes');
  tmpsd2 = ft_checkdata(ft_selectdata(source, 'rpt', indx1b), 'sourcerepresentation', 'new', 'haspow', 'yes');
  tmpsd3 = ft_checkdata(ft_selectdata(source, 'rpt', indx2a), 'sourcerepresentation', 'new', 'haspow', 'yes');
  tmpsd4 = ft_checkdata(ft_selectdata(source, 'rpt', indx2b), 'sourcerepresentation', 'new', 'haspow', 'yes');
 
  tmpsd1 = rmfield(tmpsd1, 'mom');
  tmpsd1 = rmfield(tmpsd1, 'momdimord');
  tmpsd2 = rmfield(tmpsd2, 'mom');
  tmpsd2 = rmfield(tmpsd2, 'momdimord');
  tmpsd3 = rmfield(tmpsd3, 'mom');
  tmpsd3 = rmfield(tmpsd3, 'momdimord');
  tmpsd4 = rmfield(tmpsd4, 'mom');
  tmpsd4 = rmfield(tmpsd4, 'momdimord');
  
  n1 = numel(tmpsd1.cumtapcnt);
  n2 = numel(tmpsd2.cumtapcnt);
  n3 = numel(tmpsd3.cumtapcnt);
  n4 = numel(tmpsd4.cumtapcnt);
  
  tmpsd1.pow = trimmean(log10(tmpsd1.pow),0.2,2);
  tmpsd1.powdimord = 'pos';
  tmpsd2.pow = trimmean(log10(tmpsd2.pow),0.2,2);
  tmpsd2.powdimord = 'pos';
  tmpsd3.pow = trimmean(log10(tmpsd3.pow),0.2,2);
  tmpsd3.powdimord = 'pos';
  tmpsd4.pow = trimmean(log10(tmpsd4.pow),0.2,2);
  tmpsd4.powdimord = 'pos';

  stat.stat1 = tmpsd1.pow;
  stat.stat2 = tmpsd2.pow;
  stat.stat3 = tmpsd3.pow;
  stat.stat4 = tmpsd4.pow;

  stat.inside  = sdorig.inside;
  stat.outside = sdorig.outside;
  stat.dim     = dim;
  stat.fwhm    = fwhm(:);
  stat.df1     = n1;
  stat.df2     = n2;
  stat.df3     = n3;
  stat.df4     = n4;

  if k==1,
    statl = stat;
  elseif k==2,
    statr = stat;
  end
end
