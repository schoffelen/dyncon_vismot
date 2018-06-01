function doSourcestatisticsPrePst(band,doplot,nrand)

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
for k = 1:length(d)
  if ~isempty(strfind(d(k).name,band)) && isempty(strfind(d(k).name,'grandavg')),
    fname = d(k).name;
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
    names{end+1} = fname(1:5);
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
stat.stat2(   nsubj+1:2*nsubj,:,:) = 0;
stat13.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat42.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat13b.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat42b.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat13x.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat42x.stat2(nsubj+1:2*nsubj,:,:) = 0;
clear s s13 s42 s13b s42b s13x s42x;

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'pooledT';
cfg.parameter = 'stat2';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.frequency = [stat.freq(1) stat.freq(end)];
cfg.avgoverfreq = 'yes';
st          = sourcestatistics(cfg, stat);
st13        = sourcestatistics(cfg, stat13);
st42        = sourcestatistics(cfg, stat42);
st13b       = sourcestatistics(cfg, stat13b);
st42b       = sourcestatistics(cfg, stat42b);
st13x       = sourcestatistics(cfg, stat13x);
st42x       = sourcestatistics(cfg, stat42x);

dim     = pos2dim3d(st.pos);
inside  = st.inside;
outside = st.outside;
if doplot,
  cd /analyse/4/Project0030/figures/sourcedata/4Dprepst
  for k = 1:length(names)
    fname = names{k};
    tmp   = reshape(mean(stat.stat2(k,:,:),3),dim);
    tmp(outside) = nan;
    tmp(1) = -max(abs(tmp(:)));
    tmp(2) =  max(abs(tmp(:)));
    figure;volplotJM(tmp,'montage');
    cd prepst
    print(gcf,'-dpng',[fname,band,'PrePst']);
    cd ..
    close
    tmp1 = reshape(mean(stat42.stat2(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13.stat2(k,:,:),3),dim);
    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
    tmp1(outside) = nan;
    tmp1(1) = -max(abs(tmp1(:)));
    tmp1(2) =  max(abs(tmp1(:)));
    tmp2(outside) = nan;
    tmp2(1) = -max(abs(tmp2(:)));
    tmp2(2) =  max(abs(tmp2(:)));
    tmp(outside) = nan;
    tmp(1)  = -max(abs(tmp(:)));
    tmp(2)  =  max(abs(tmp(:)));
    figure;
    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
    cd congruency
    print(gcf,'-dpng',[fname,band,'Congruency']);
    cd ..
    close
    tmp1 = reshape(mean(stat42x.stat2(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13x.stat2(k,:,:),3),dim);
    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
    tmp1(outside) = nan;
    tmp1(1) = -max(abs(tmp1(:)));
    tmp1(2) =  max(abs(tmp1(:)));
    tmp2(outside) = nan;
    tmp2(1) = -max(abs(tmp2(:)));
    tmp2(2) =  max(abs(tmp2(:)));
    tmp(outside) = nan;
    tmp(1)  = -max(abs(tmp(:)));
    tmp(2)  =  max(abs(tmp(:)));
    figure;
    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
    cd congruencyRelative
    print(gcf,'-dpng',[fname,band,'Congruency']);
    cd ..
    close
    tmp1 = reshape(mean(stat42b.stat2(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13b.stat2(k,:,:),3),dim);
    tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
    tmp1(outside) = nan;
    tmp1(1) = -max(abs(tmp1(:)));
    tmp1(2) =  max(abs(tmp1(:)));
    tmp2(outside) = nan;
    tmp2(1) = -max(abs(tmp2(:)));
    tmp2(2) =  max(abs(tmp2(:)));
    tmp(outside) = nan;
    tmp(1)  = -max(abs(tmp(:)));
    tmp(2)  =  max(abs(tmp(:)));
    figure;
    subplot('position',[0.26 0.01 0.48 0.48]);volplotJM(tmp,'montage'); colorbar off
    subplot('position',[0.01 0.51 0.48 0.48]);volplotJM(tmp2,'montage'); colorbar off
    subplot('position',[0.51 0.51 0.48 0.48]);volplotJM(tmp1,'montage'); colorbar off
    cd congruencyBaseline
    print(gcf,'-dpng',[fname,band,'Congruency']);
    cd ..
    close
  end
  
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
  print(gcf,'-dpng',[band,'AvgCongruency']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyRight']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyLeft']);
  close
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
  print(gcf,'-dpng',[band,'AvgCongruencyRelative']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyRightRelative']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyLeftRelative']);
  close
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
  print(gcf,'-dpng',[band,'AvgCongruencyBaseline']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyRightBaseline']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgCongruencyLeftBaseline']);
  close
  tmp = reshape(st.stat,dim);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'AvgPrePst']);
  close
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
  dum = zeros(dim);
  for k = 1:length(stat42.freq)
    for m = 1:nsubj
      dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2(m,:,k) = dum(:);
      dum(:) = (reshape(stat42b.stat2(m,:,k),dim)+flipdim(reshape(stat13b.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s2.stat2(m,:,k) = dum(:);
    end
  end
  sel = find(~ismember(names,{'GAR12' 'BKA01'}));
  sel = [sel sel+nsubj];
  cfg.design = cfg.design(:,sel);
  s   = selectdata(s,'rpt',sel);
  s2  = selectdata(s2,'rpt',sel);
  cfg.numrandomization = nrand;
  cfg.correctm         = 'max';
  statC = sourcestatistics(cfg, s);
  statB = sourcestatistics(cfg, s2);
  save(['grandavg_',band],'statC','statB','st');
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
