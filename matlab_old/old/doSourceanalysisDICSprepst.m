function [stat, stat13, stat42, stat13b, stat42b, stat13x, stat42x] = doSourceanalysisDICSprepst(subject, frequency, flag)

%do sourceanalysis using all baseline and post stimulus data (number of samples for each condition's baseline and post
%stimulus interval are matched, and number of samples across conditions which are contrasted are matched)
%spatial filters are computed on all trials concatenated
%output:
%stat: act vs baseline depsamplesT
%stat13 stat42: act congruency effect
%stat13b stat42b: baselines congruency effect
%stat13x stat42x: relative change (log10(act)-log10(bas)) congruency effect

if nargin<3,
  flag = 0;
end

cd([subject.pathname,'freq/']);
if flag==0 && frequency<37,
  load([subject.name,'mtmfft004b']);
  %load([subject.name,'hanning'], 'freqpst');
elseif flag==0,
  load([subject.name,'mtmfft012b']);
elseif flag==1 && frequency<37,
  load([subject.name,'mtmfft_aligned004b']);
elseif flag==2 && frequency<37,
  load([subject.name,'mtmfft_aligned004bc']);
elseif flag==1,
  load([subject.name,'mtmfft_aligned012b']);
elseif flag==2,
  load([subject.name,'mtmfft_aligned012bc']);
end

%concatenate data and keep track of original trial numbers
gradlabel = freqpst{1}.grad.label;
ok = zeros(248,1);
for k = 1:5
  [a,b] = match_str(gradlabel, freqpst{k}.label);
  ok(a) = ok(a)+1;
  ntrl(1,k) = length(freqpre{k}.cumtapcnt);
  ntrl(1,k+5) = length(freqpst{k}.cumtapcnt);
end
oklabel = gradlabel(ok==5);
for k = 1:5
  warning off
  freqpst{k} = struct2double(freqpst{k});
  freqpre{k} = struct2double(freqpre{k});
  warning on;
  freqpst{k} = selectdata(freqpst{k}, 'channel', oklabel);
  freqpre{k} = selectdata(freqpre{k}, 'channel', oklabel);
end

freq{1} = selectdata(freqpre{1:5}, 'param', 'fourierspctrm');
freq{2} = selectdata(freqpst{1:5}, 'param', 'fourierspctrm');
freq    = selectdata(freq{:}, 'param', 'fourierspctrm');
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
%1:5 are the baselines for conditions 1:5
%6:10 are the post stim intervals for conditions 1:5

freqx = selectdata(freq, 'foilim', frequency);

%project the leadfields onto the optimal orientation
cfg.grid.filter = source.avg.filter;
for k = 1:length(source.inside)
  indx = source.inside(k);
  cfg.grid.leadfield{indx} = newgrid.leadfield{indx}*source.avg.ori{indx};
end
cfg.method      = 'pcc';
cfg.keepmom     = 'yes';

%project the single trials through the filters
source          = sourceanalysis(cfg, freqx);
addpath /home/jan/projects/ccc/3D
sd              = estimate_fwhm_tetra(source);
sd              = checkdata(sd, 'hasdimord', 'yes', 'sourcedimord', 'rpt_pos');
clear source;

fwhm          = sd.fwhm;
krn           = compute_kernel(sd, 'truncate', 5e-5);

sd = rmfield(sd, 'fwhm');
sd = rmfield(sd, 'leadfield');
sd = rmfield(sd, 'rough');
sd = rmfield(sd, 'q');

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'indepsamplesT';
cfg.parameter = 'pow';

%activation versus baseline on all trials
cfg.statistic = 'depsamplesT';
cfg.design    = ones(1,size(sd.pow,1));
cfg.design(1:etrl(5)) = 2;
cfg.design(2,:) = [1:etrl(5) 1:etrl(5)];
cfg.ivar      = 1;
cfg.uvar      = 2;
stat          = sourcestatistics(cfg, sd);
dim           = pos2dim3d(stat.pos);
inside        = stat.inside;
stat.stat2    = smooth_vol(stat.stat,krn,dim,inside);

%activation condition 1 vs 3
cfg.statistic = 'indepsamplesT';
cfg           = rmfield(cfg,'uvar');
tmpsd = selectdata(sd, 'rpt', [btrl(6):etrl(6) btrl(8):etrl(8)]);
n1    = numel(btrl(6):etrl(6));
n2    = numel(btrl(8):etrl(8));
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat13     = sourcestatistics(cfg, tmpsd);
stat13.stat2 = smooth_vol(stat13.stat,krn,dim,inside);
stat13.cfg   = [];

%activation condition 4 vs 2
tmpsd = selectdata(sd, 'rpt', [btrl(9):etrl(9) btrl(7):etrl(7)]);
n1    = numel(btrl(9):etrl(9));
n2    = numel(btrl(7):etrl(7));
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat42     = sourcestatistics(cfg, tmpsd);
stat42.stat2 = smooth_vol(stat42.stat,krn,dim,inside);
stat42.cfg   = [];

%bsl condition 1 vs 3
tmpsd = selectdata(sd, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
n1    = numel(btrl(1):etrl(1));
n2    = numel(btrl(3):etrl(3));
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat13b    = sourcestatistics(cfg, tmpsd);
stat13b.stat2 = smooth_vol(stat13b.stat,krn,dim,inside);
stat13b.cfg   = [];

%bsl condition 4 vs 2
tmpsd = selectdata(sd, 'rpt', [btrl(4):etrl(4) btrl(2):etrl(2)]);
n1    = numel(btrl(4):etrl(4));
n2    = numel(btrl(2):etrl(2));
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat42b    = sourcestatistics(cfg, tmpsd);
stat42b.stat2 = smooth_vol(stat42b.stat,krn,dim,inside);
stat42b.cfg   = [];

%relative change condition 1 vs 3
tmpsd  = selectdata(sd, 'rpt', [btrl(6):etrl(6) btrl(8):etrl(8)]);
tmpsdb = selectdata(sd, 'rpt', [btrl(1):etrl(1) btrl(3):etrl(3)]);
n1     = numel(btrl(6):etrl(6));
n2     = numel(btrl(8):etrl(8));
tmpsd.pow  = log10(tmpsd.pow)-log10(tmpsdb.pow);
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat13x = sourcestatistics(cfg, tmpsd);
stat13x.stat2 = smooth_vol(stat13x.stat,krn,dim,inside);
stat13x.cfg   = [];

%relative change condition 4 vs 2
tmpsd  = selectdata(sd, 'rpt', [btrl(9):etrl(9) btrl(7):etrl(7)]);
tmpsdb = selectdata(sd, 'rpt', [btrl(4):etrl(4) btrl(2):etrl(2)]);
n1     = numel(btrl(9):etrl(9));
n2     = numel(btrl(7):etrl(7));
tmpsd.pow  = log10(tmpsd.pow)-log10(tmpsdb.pow);
cfg.design = [ones(1,n1) ones(1,n2)*2];
stat42x = sourcestatistics(cfg, tmpsd);
stat42x.stat2 = smooth_vol(stat42x.stat,krn,dim,inside);
stat42x.cfg   = [];
