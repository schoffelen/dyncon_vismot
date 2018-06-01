d('/analyse/4/Project0030/source/4Dglm');

d = dir;
d = d(3:end);
cnt = 0;

load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = 1:length(d)
  fname = d(k).name;
  fprintf('loading %s\n', fname);
  load(fname);
  stat13.pos = grid.pos;
  stat42.pos = grid.pos;
  cnt      = cnt+1;
  s13{cnt} = stat13;
  s42{cnt} = stat42;
end
stat13 = selectdata(s13{:}, 'param', 'stat2');
stat42 = selectdata(s42{:}, 'param', 'stat2');
nsubj  = size(stat13.stat2,1);
%stat13.stat2(nsubj+1:2*nsubj,:,:) = 0;
%stat42.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat13.stat2(nsubj+1:2*nsubj,:,:) = repmat(nanmean(stat13.stat2(1:nsubj,:,:),2),[1 33480 1]);
stat42.stat2(nsubj+1:2*nsubj,:,:) = repmat(nanmean(stat42.stat2(1:nsubj,:,:),2),[1 33480 1]);

stat13.dim(1) = 2*nsubj;
stat42.dim(1) = 2*nsubj;

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'pooledT';
%cfg.statistic = 'depsamplesT';
cfg.parameter = 'stat2';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.uvar      = 2;
s13           = sourcestatistics(cfg, stat13);
s42           = sourcestatistics(cfg, stat42);

inside = logical(zeros(size(s13.stat)));
inside(s13.inside,:) = 1;
s13.inside = inside;
s42.inside = inside;


s13.dim = s13.dim(2:end);
s42.dim = s42.dim(2:end);
s13.dimord = 'pos_freq';
s42.dimord = 'pos_freq';

%cfgp = [];
%cfgp.method       = 'ortho';
%cfgp.interactive  = 'yes';
%cfgp.funparameter = 'stat';
%figure;sourceplot(cfgp, s13);

addpath /home/jan/matlab/spm2;
s   = stat42;
dim = pos2dim3d(stat42.pos);
dum = zeros(dim);
for k = 1:length(stat42.freq)
  for m = 1:nsubj
    dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
    spm_smooth(dum,dum,1.5);
    s.stat2(m,:,k) = dum(:);
  end
end
cfg.numrandomization = 2000;
cfg.correctm         = 'max';
stat = sourcestatistics(cfg, s);
stat.dimord = 'pos_freq';
%stat.dim    = stat.dim(2:end);

tmpinside   = logical(zeros(size(stat.stat)));
tmpinside(stat.inside,:) = 1;
stat.inside = tmpinside;

%mri = read_mri('/home/jan/matlab/mri/templateMRI.mnc');
%cfgi.parameter = 'stat';
%i1  = sourceinterpolate(cfgi, stat, mri);
