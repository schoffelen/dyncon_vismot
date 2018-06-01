function doSourcestatisticsBehav1(band,doplot,nrand)

subjinfo;

if nargin<3,
  nrand = 0;
end

if nargin<2,
  doplot = 0;
end

cd('/analyse/4/Project0030/source/4Dprepst');

d = dir;
d = d(3:end);
cnt = 0;
names = {};
load('/home/jan/matlab/mri/templategrid6mm.mat');
for k = [1:3 5:6 8:16 18:20]
  %Exclude GAR12 (weird behaviour) and BKA01 (bad data)
  x1  = cellfun(@isempty,strfind({d(:).name},band))==0;
  x2  = cellfun(@isempty,strfind({d(:).name},SUBJ(k).name))==0;
  sel = find(x1 & x2);
  
  if ~isempty(sel),
  fname = d(sel).name;
  fprintf('loading %s\n', fname);
  load(fname);
  stat.pos    = grid.pos;
  stat13.pos  = grid.pos;
  stat42.pos  = grid.pos;
  stat13b.pos = grid.pos;
  stat42b.pos = grid.pos;
  stat13x.pos = grid.pos;
  stat42x.pos = grid.pos;
  cnt       = cnt+1;
  s{cnt}    = stat;
  s13{cnt}  = stat13;
  s42{cnt}  = stat42;
  s13b{cnt} = stat13b;
  s42b{cnt} = stat42b;
  s13x{cnt} = stat13x;
  s42x{cnt} = stat42x;
  rt(cnt,:) = SUBJ(k).rttrim;
  names{cnt} = SUBJ(k).name;
  end
end
stat    = selectdata(s{:},    'param', 'stat2');
stat13  = selectdata(s13{:},  'param', 'stat2');
stat42  = selectdata(s42{:},  'param', 'stat2');
stat13b = selectdata(s13b{:}, 'param', 'stat2');
stat42b = selectdata(s42b{:}, 'param', 'stat2');
stat13x = selectdata(s13x{:}, 'param', 'stat2');
stat42x = selectdata(s42x{:}, 'param', 'stat2');
nsubj   = size(stat.stat2,1);
clear s s13 s42 s13b s42b s13x s42x;

stat   = selectdata(stat, 'avgoverfreq', 'yes');
stat13 = selectdata(stat13, 'avgoverfreq', 'yes');
stat13b = selectdata(stat13b, 'avgoverfreq', 'yes');
stat13x = selectdata(stat13x, 'avgoverfreq', 'yes');
stat42 = selectdata(stat42, 'avgoverfreq', 'yes');
stat42b = selectdata(stat42b, 'avgoverfreq', 'yes');
stat42x = selectdata(stat42x, 'avgoverfreq', 'yes');

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'spearman';
cfg.precondition = 'before';
cfg.parameter = 'stat2';
cfg.design    = [rt(:,1)-rt(:,3)]';
cfg.ivar      = 1;
%st          = sourcestatistics(cfg, stat);
st13        = sourcestatistics(cfg, stat13);
st13b       = sourcestatistics(cfg, stat13b);
st13x       = sourcestatistics(cfg, stat13x);
cfg.design  = [rt(:,4)-rt(:,2)]';
st42        = sourcestatistics(cfg, stat42);
st42b       = sourcestatistics(cfg, stat42b);
st42x       = sourcestatistics(cfg, stat42x);

dim     = pos2dim3d(st13.pos);
inside  = st13.inside;
outside = st13.outside;
insidevol = zeros(dim);
insidevol(inside) = 1;
inside = find(insidevol+flipdim(insidevol,1)==2);
outside = setdiff(1:prod(dim),inside);
statall = stat13;
flipop  = reshape(1:prod(dim),dim);
flipop  = flipdim(flipop,1);
flipop  = flipop(:);
statall.stat2 = double([stat42.stat2;stat13.stat2(:,flipop)]);
statall.inside = inside;
statall.outside = outside;
statall.dim(1)= 2*nsubj;
cfg.design = [rt(:,4)-rt(:,2);rt(:,1)-rt(:,3)]';
nrand2 = 0;
cfg.numrandomization = nrand2;
stall   = sourcestatistics(cfg,statall);
statall.stat2 = double([stat42b.stat2;stat13b.stat2(:,flipop)]);
stallb  = sourcestatistics(cfg,statall);
statall.stat2 = double([stat42x.stat2;stat13x.stat2(:,flipop)]);
stallx  = sourcestatistics(cfg,statall);

if doplot,
  cd /analyse/4/Project0030/figures/sourcedata/prepstBehav1
  
  %plot the averages
  tmp1 = reshape(st42.stat,dim);
  tmp2 = reshape(st13.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'Behav1pooled']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'Behav1Right']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'Behav1Left']);
  close
  tmp  = reshape(stall.stat,dim);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'pooledBehav1']);
  close
  if nrand2>0,
    tmp = reshape(-log10(stall.stat+1/(10*nrand)),dim);
    tmp(1)=-0.5;
    figure;volplotJM(tmp,'montage');
    print(gcf,'-dpng',[band,'pooledBehav1P']);
    close
  end
  tmp1 = reshape(st42x.stat,dim);
  tmp2 = reshape(st13x.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'Behav1pooledRelative']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'Behav1RightRelative']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'Behav1LeftRelative']);
  close
  tmp  = reshape(stallx.stat,dim);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'pooledBehav1Relative']);
  close
  if nrand2>0,
    tmp = reshape(-log10(stallx.stat+1/(10*nrand)),dim);
    tmp(1)=-0.5;
    figure;volplotJM(tmp,'montage');
    print(gcf,'-dpng',[band,'pooledBehav1RelativeP']);
    close
  end
  tmp1 = reshape(st42b.stat,dim);
  tmp2 = reshape(st13b.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'Behav1pooledBaseline']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'Behav1RightBaseline']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'Behav1LeftBaseline']);
  close
  tmp  = reshape(stallb.stat,dim);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'pooledBehav1Baseline']);
  close
  if nrand2>0,
    tmp = reshape(-log10(stallb.stat+1/(10*nrand)),dim);
    tmp(1)=-0.5;
    figure;volplotJM(tmp,'montage');
    print(gcf,'-dpng',[band,'pooledBehav1BaselineP']);
    close
  end
end

%s13.dim = s13.dim(2:end);
%s42.dim = s42.dim(2:end);
%s13.dimord = 'pos_freq';
%s42.dimord = 'pos_freq';
%
%%cfgp = [];
%%cfgp.method       = 'ortho';
%%cfgp.interactive  = 'yes';
%%cfgp.funparameter = 'stat';
%%figure;sourceplot(cfgp, s13);
%

if nrand>0,
  addpath /home/jan/matlab/spm2;
  s   = stat42;
  s2  = stat42b;
  s3  = stat42x;
  dum = zeros(dim);
  for k = 1:length(stat42.freq)
    for m = 1:nsubj
      dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2(m,:,k) = dum(:);
      dum(:) = (reshape(stat42b.stat2(m,:,k),dim)+flipdim(reshape(stat13b.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s2.stat2(m,:,k) = dum(:);
      dum(:) = (reshape(stat42x.stat2(m,:,k),dim)+flipdim(reshape(stat13x.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s3.stat2(m,:,k) = dum(:);
    end
  end
  sel = find(~ismember(names,{'GAR12' 'BKA01'}));
  %sel = [sel sel+nsubj];
  cfg.design = cfg.design(:,sel);
  s   = selectdata(s,'rpt',sel);
  s2  = selectdata(s2,'rpt',sel);
  s3  = selectdata(s3,'rpt',sel);
  cfg.numrandomization = nrand;
  %cfg.correctm         = 'max';
  statC = sourcestatistics(cfg, s);
  statB = sourcestatistics(cfg, s2);
  statX = sourcestatistics(cfg, s3);
  cd /analyse/4/Project0030/source/4DprepstBehav1
  save(['grandavg_',band,'Behav1'],'statC','statB','statX');
end

%stat.dimord = 'pos_freq';
%%stat.dim    = stat.dim(2:end);
%
%tmpinside   = logical(zeros(size(stat.stat)));
%tmpinside(stat.inside,:) = 1;
%stat.inside = tmpinside;
%
%%mri = read_mri('/home/jan/matlab/mri/templateMRI.mnc');
%%cfgi.parameter = 'stat';
%%i1  = sourceinterpolate(cfgi, stat, mri);
