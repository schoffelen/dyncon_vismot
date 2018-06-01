function doSourcestatisticsPrePstGLM(band,doplot,nrand)

if nargin<3,
  nrand = 0;
end

if nargin<2,
  doplot = 0;
end

cd('/analyse/4/Project0030/source/4DprepstGLM');

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
    stat13.pos  = grid.pos;
    stat42.pos  = grid.pos;
    stat13b.pos = grid.pos;
    stat42b.pos = grid.pos;
    cnt       = cnt+1;
    s13{cnt}  = stat13;
    s42{cnt}  = stat42;
    s13b{cnt} = stat13b;
    s42b{cnt} = stat42b;
    names{end+1} = fname(1:5);
  end
end
stat13  = selectdata(s13{:},  'param', {'stat2' 'stat2rt'});
stat42  = selectdata(s42{:},  'param', {'stat2' 'stat2rt'});
stat13b = selectdata(s13b{:}, 'param', {'stat2' 'stat2rt'});
stat42b = selectdata(s42b{:}, 'param', {'stat2' 'stat2rt'});
nsubj   = size(stat13.stat2,1);
stat13.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat42.stat2( nsubj+1:2*nsubj,:,:) = 0;
stat13b.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat42b.stat2(nsubj+1:2*nsubj,:,:) = 0;
stat13.stat2rt( nsubj+1:2*nsubj,:,:) = 0;
stat42.stat2rt( nsubj+1:2*nsubj,:,:) = 0;
stat13b.stat2rt(nsubj+1:2*nsubj,:,:) = 0;
stat42b.stat2rt(nsubj+1:2*nsubj,:,:) = 0;
clear s13 s42 s13b s42b;

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic = 'pooledT';
cfg.parameter = 'stat2';
cfg.design    = [ones(1,nsubj) ones(1,nsubj)*2;1:nsubj 1:nsubj];
cfg.ivar      = 1;
cfg.frequency = [stat13.freq(1) stat13.freq(end)];
cfg.avgoverfreq = 'yes';
st13        = sourcestatistics(cfg, stat13);
st42        = sourcestatistics(cfg, stat42);
st13b       = sourcestatistics(cfg, stat13b);
st42b       = sourcestatistics(cfg, stat42b);
cfg.parameter = 'stat2rt';
st13rt        = sourcestatistics(cfg, stat13);
st42rt        = sourcestatistics(cfg, stat42);
st13brt       = sourcestatistics(cfg, stat13b);
st42brt       = sourcestatistics(cfg, stat42b);

dim     = pos2dim3d(st13.pos);
inside  = st13.inside;
outside = st13.outside;
if doplot,
  cd /analyse/4/Project0030/figures/sourcedata/4DprepstGLM
  for k = 1:length(names)
    fname = names{k};
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
    tmp1 = reshape(mean(stat42.stat2rt(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13.stat2rt(k,:,:),3),dim);
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
    cd reactiontime
    print(gcf,'-dpng',[fname,band,'Reactiontime']);
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
    tmp1 = reshape(mean(stat42b.stat2rt(k,:,:),3),dim);
    tmp2 = reshape(mean(stat13b.stat2rt(k,:,:),3),dim);
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
    cd reactiontimeBaseline
    print(gcf,'-dpng',[fname,band,'Reactiontime']);
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
  tmp1 = reshape(st42rt.stat,dim);
  tmp2 = reshape(st13rt.stat,dim);
  tmp  = (tmp1+flipdim(tmp2,1))./sqrt(2);
  tmp(1) = -max(abs(tmp(:)));
  tmp(2) =  max(abs(tmp(:)));
  tmp1(1) = -max(abs(tmp1(:)));
  tmp1(2) =  max(abs(tmp1(:)));
  tmp2(1) = -max(abs(tmp2(:)));
  tmp2(2) =  max(abs(tmp2(:)));
  figure;volplotJM(tmp,'montage');
  print(gcf,'-dpng',[band,'AvgReactiontime']);
  close
  figure;volplotJM(tmp1,'montage');
  print(gcf,'-dpng',[band,'AvgReactiontimeRight']);
  close
  figure;volplotJM(tmp2,'montage');
  print(gcf,'-dpng',[band,'AvgReactiontimeLeft']);
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
  tmp1 = reshape(st42brt.stat,dim);
  tmp2 = reshape(st13brt.stat,dim);
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
end

if nrand>0,
  addpath /home/jan/matlab/spm2;
  s   = stat42;
  s2  = stat42b;
  s.stat2l = s.stat2;
  s.stat2r = s.stat2;
  s2.stat2l = s2.stat2;
  s2.stat2r = s2.stat2;
  dum = zeros(dim);
  for k = 1:length(stat42.freq)
    for m = 1:nsubj
      dum(:) = reshape(stat42.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2r(m,:,k) = dum(:);
      dum(:) = reshape(stat13.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s.stat2l(m,:,k) = dum(:);
      dum(:) = (reshape(stat42.stat2(m,:,k),dim)+flipdim(reshape(stat13.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2(m,:,k) = dum(:);
      dum(:) = (reshape(stat42b.stat2(m,:,k),dim)+flipdim(reshape(stat13b.stat2(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s2.stat2(m,:,k) = dum(:);
      dum(:) = reshape(stat42b.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s2.stat2r(m,:,k) = dum(:);
      dum(:) = reshape(stat13b.stat2(m,:,k),dim);
      spm_smooth(dum,dum,1.5);
      s2.stat2l(m,:,k) = dum(:);
      dum(:) = (reshape(stat42.stat2rt(m,:,k),dim)+flipdim(reshape(stat13.stat2rt(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s.stat2rt(m,:,k) = dum(:);
      dum(:) = (reshape(stat42b.stat2rt(m,:,k),dim)+flipdim(reshape(stat13b.stat2rt(m,:,k),dim),1))./sqrt(2);
      spm_smooth(dum,dum,1.5);
      s2.stat2rt(m,:,k) = dum(:);
    end
  end
  sel = find(~ismember(names,{'GAR12' 'BKA01'}));
  sel = [sel sel+nsubj];
  cfg.design = cfg.design(:,sel);
  s   = selectdata(s,'rpt',sel);
  s2  = selectdata(s2,'rpt',sel);
  cfg.numrandomization = nrand;
  cfg.correctm         = 'max';
  cfg.parameter        = 'stat2';
  statC = sourcestatistics(cfg, s);
  statB = sourcestatistics(cfg, s2);
  cfg.parameter = 'stat2l';
  statCl = sourcestatistics(cfg, s);
  statBl = sourcestatistics(cfg, s2);
  cfg.parameter = 'stat2r';
  statCr = sourcestatistics(cfg, s);
  statBr = sourcestatistics(cfg, s2);
  cd /analyse/4/Project0030/source/4DprepstGLM
  save(['grandavg_',band],'statC','statB','statCl','statCr','statBl','statBr');
  cfg.parameter        ='stat2rt';
  statC = sourcestatistics(cfg, s);
  statB = sourcestatistics(cfg, s2);
  save(['grandavgRT_',band],'statC','statB');
end
