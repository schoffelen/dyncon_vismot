function [stat13, stat42] = doSourceanalysisDICSglm(subject, frequency, flag)

if nargin<3,
  flag = 0;
end

cd([subject.pathname,'freq/']);
if flag==0 && frequency<37,
  load([subject.name,'mtmfft004'], 'freqpst');
  %load([subject.name,'hanning'], 'freqpst');
elseif flag==0,
  load([subject.name,'mtmfft012'], 'freqpst');
elseif flag==1 && frequency<37,
  load([subject.name,'mtmfft_aligned004'], 'freqpst');
elseif flag==1,
  load([subject.name,'mtmfft_aligned012'], 'freqpst');
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
for k = 1:4
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(k) = length(freqpst{k}.cumtapcnt);
end
oklabel = gradlabel(ok==4);
oktrl   = {randperm(ntrl(1)) randperm(ntrl(2)) randperm(ntrl(3)) randperm(ntrl(4))};
mintrl  = min(ntrl);
for k = 1:4
  oktrl{k}   = oktrl{k}(1:mintrl);
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel, 'rpt', oktrl{k});
end
ntrl(:) = mintrl;

freq = selectdata(freqpst{1:4}, 'param', 'fourierspctrm');
clear freqpst

%load grid and prune leadfields
[a,b] = match_str(freq.label, gradlabel);
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
cfg.lambda     = '15%';
cfg.feedback   = 'none';
%[u,s,v] = svd(tlckcovpre.cov);
%cfg.subspace   = pinv(sqrtm(tlckcovpre.cov+eye(length(tlckcovpre.label))*0.01*s(1,1)));
%cfg.lambda     = 1;
source         = sourceanalysis(cfg, freq);

strl = cumsum([0 ntrl]);
btrl = strl+1;
etrl = strl(2:end);

freq13 = selectdata(freq, 'foilim', frequency);
freq13 = selectdata(freq13, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
freq24 = selectdata(freq, 'foilim', frequency);
freq24 = selectdata(freq24, 'rpt', [btrl(2):etrl(2) btrl(4):etrl(4)]);
%clear freq

cfg.grid.filter = source.avg.filter;
for k = 1:length(source.inside)
  indx = source.inside(k);
  cfg.grid.leadfield{indx} = newgrid.leadfield{indx}*source.avg.ori{indx};
end
cfg.method      = 'pcc';
cfg.keepmom     = 'yes';
source13        = sourceanalysis(cfg, freq13);
source24        = sourceanalysis(cfg, freq24);
cumsumcnt13     = freq13.cumsumcnt;
cumsumcnt24     = freq24.cumsumcnt;
clear source freq13 freq24;

addpath /home/jan/projects/ccc/3D
sd13 = estimate_fwhm_tetra(source13);
sd24 = estimate_fwhm_tetra(source24);

sd13 = checkdata(sd13, 'hasdimord','yes','sourcedimord','rpt_pos');
sd24 = checkdata(sd24, 'hasdimord','yes','sourcedimord','rpt_pos');

%cfg = [];
%cfg.keeptrials = 'yes';
%cfg.keepmom    = 'no';
%sd13 = sourcedescriptives(cfg, sd13);
%sd24 = sourcedescriptives(cfg, sd24);

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'glmJM';
cfg.glm.statistic = 'T';
cfg.glm.contrast  = [1 0];
cfg.parameter = 'pow';
cfg.design    = [ones(1,ntrl(1)) ones(1,ntrl(3))*2;cumsumcnt13(:)'];
stat13        = sourcestatistics(cfg, sd13);
krn           = compute_kernel(sd13, 'truncate', 5e-5);
stat13.stat2  = smooth_vol(stat13.stat,krn,pos2dim3d(stat13.pos),stat13.inside);
clear krn;

cfg.design    = [ones(1,ntrl(2))*2 ones(1,ntrl(4));cumsumcnt24(:)'];
stat42        = sourcestatistics(cfg, sd24);
krn           = compute_kernel(sd24, 'truncate', 5e-5);
stat42.stat2  = smooth_vol(stat42.stat,krn,pos2dim3d(stat42.pos),stat42.inside);
