subjinfo;
load /home/jan/projects/visuomotor/matlab/lrplist
cd /analyse/4/Project0030/freq
cnt = 0;
for k = 1:length(SUBJ)
  %fname = [SUBJ(k).name,'lrp'];
  fname = [SUBJ(k).name,'planar004'];
  try,
    load(fname, 'stat');
    fprintf('processing %s\n', fname);
    cnt = cnt+1;
  end
  [a1,b1]   = match_str(channelcmb(:,1),stat(3).label);
  [a2,b2]   = match_str(channelcmb(:,2),stat(2).label);
  [int,ia1,ia2] = intersect(a1,a2);
  tmp       = stat(3);
  tmp.label = stat(3).label(b1(ia1));
  tmp.stat  = (stat(3).stat(b1(ia1),:) + stat(2).stat(b2(ia2),:))./sqrt(2);
  dat{cnt}  = tmp;
  [a1,b1]   = match_str(channelcmb(:,1),stat(5).label);
  [a2,b2]   = match_str(channelcmb(:,2),stat(4).label);
  [int,ia1,ia2] = intersect(a1,a2);
  tmp       = stat(5);
  tmp.label = stat(5).label(b1(ia1));
  tmp.stat  = (stat(5).stat(b1(ia1),:) + stat(4).stat(b2(ia2),:))./sqrt(2);
  datb{cnt}  = tmp;
end

cfg.layout = '4D248.lay';
lay  = prepare_layout(cfg);
cnt  = zeros(248,1);
cntb = zeros(248,1);
for k = 1:length(dat)
  [a,b] = match_str(lay.label, dat{k}.label);
  cnt(a) = cnt(a)+1;
  [a,b] = match_str(lay.label, datb{k}.label);
  cntb(a) = cntb(a)+1;
end

sel     = find(cnt==length(dat));
selb    = find(cntb==length(datb));
channel = lay.label(sel);
for k = 1:length(dat)
  dat{k}.powspctrm  = dat{k}.stat;
  datb{k}.powspctrm = datb{k}.stat;
end
for k = 1:length(dat)
  dat2{k}  = selectdata(dat{k},'channel',channel);
  dat2b{k} = selectdata(datb{k},'channel',channel);
end
dat    = dat2;clear dat2;
datb   = dat2b; clear dat2b;
for k = 1:length(dat)
  dat{k} = rmfield(dat{k}, 'stat');
  dat{k} = rmfield(dat{k}, 'prob');
  dat{k} = rmfield(dat{k}, 'mask');
  datb{k} = rmfield(datb{k}, 'stat');
  datb{k} = rmfield(datb{k}, 'prob');
  datb{k} = rmfield(datb{k}, 'mask');
end
alldat  = selectdata(dat{:},'param','powspctrm');
alldatb = selectdata(datb{:},'param','powspctrm');
nsubj  = size(alldat.powspctrm,1);
alldat.powspctrm(nsubj+1:2*nsubj,:,:) = 0;
alldatb.powspctrm(nsubj+1:2*nsubj,:,:) = 0;

cd /home/jan/projects/visuomotor/matlab
cd /analyse/4/Project0030/data
load CDE04data data1;
alldat.grad          = data1.grad;
alldatb.grad         = data1.grad;
cfg                  = [];
cfg.grad             = data1.grad;
cfg.neighbourdist    = 0.037;
neighb               = neighbourselection(cfg);
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.numrandomization = 500;
cfg.statistic        = 'indepsamplesT';
%cfg.correctm         = 'cluster';
cfg.correctm         = 'no';
%cfg.clusterthreshold = 'nonparametric_common';
cfg.design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];
cfg.ivar   = 1;
cfg.uvar   = 2;
%cfg.neighbours = neighb;
stat  = freqstatistics(cfg, alldat);
statb = freqstatistics(cfg, alldatb);

cfg = [];
cfg.layout = '4D248.lay';
lay = prepare_layout(cfg);

cd /analyse/4/Project0030/figures/sensordata
cfg        = [];
cfg.layout = lay;
cfg.zparam = 'stat';
cfg.xlim   = [8 12];
cfg.zlim   = 'absmax';
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','congruencyeffectalpha');
cfg.zlim   = get(gca,'clim');
%figure;topoplotER(cfg,statb);
%
cfg.xlim   = [18 28];
cfg.zlim   = 'absmax';
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','congruencyeffectbeta');
%cfg.zlim   = get(gca,'clim');
%figure;topoplotER(cfg,statb);
%
%cfg.xlim   = [40 80];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','congruencyeffectgamma');

%cfg.xlim         = 56+[-12 12];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastgammalow12');
%cfg.xlim         = [85 100];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastgammahigh12');

%-------OLD CODE
%subjinfo;
%cnt = 0;
%for k = 1:length(SUBJ)
%  %fname = [SUBJ(k).name,'lrp'];
%  fname = [SUBJ(k).name,'lrp12'];
%  try,
%    load(fname, 'lrp');
%    cnt = cnt+1;
%  end
%  lrpx = lrp;
%  lrpx.label = lrpx.plotlabel;
%  dat{cnt} = lrpx;
%end
%
%cfg.layout = '4D248.lay';
%lay = prepare_layout(cfg);
%cnt = zeros(248,1);
%for k = 1:length(dat)
%  [a,b] = match_str(lay.label, dat{k}.label);
%  cnt(a) = cnt(a)+1;
%end
%
%sel     = find(cnt==length(dat));
%channel = lay.label(sel);
%for k = 1:length(dat)
%  dat2{k} = selectdata(dat{k},'channel',channel);
%end
%dat    = dat2;clear dat2;
%for k = 1:length(dat)
%  dat{k}.powspctrm = dat{k}.avg;
%  dat{k} = rmfield(dat{k},'avg');
%  dat{k}.dimord = 'chan_freq';
%  dat{k}.freq = dat{k}.time;
%end  
%alldat = selectdata(dat{:},'param','powspctrm');
%nsubj  = size(alldat.powspctrm,1);
%alldat.powspctrm(nsubj+1:2*nsubj,:,:) = 0;
%
%cd /home/jan/projects/visuomotor/matlab
%load leftchanlist;
%cd /analyse/4/Project0030/data
%load CDE04data data1;
%alldat.grad          = data1.grad;
%alldat               = rmfield(alldat,'time');
%cfg                  = [];
%cfg.grad             = data1.grad;
%cfg.neighbourdist    = 0.037;
%neighb               = neighbourselection(cfg);
%cfg                  = [];
%cfg.method           = 'montecarlo';
%cfg.numrandomization = 500;
%cfg.statistic        = 'pooledT';
%cfg.correctm         = 'cluster';
%cfg.clusterthreshold = 'nonparametric_common';
%cfg.design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];
%cfg.ivar   = 1;
%cfg.uvar   = 2;
%cfg.channel = channel;
%cfg.neighbours = neighb;
%stat = freqstatistics(cfg, alldat);
%
%cfg = [];
%cfg.layout = '4D248.lay';
%lay = prepare_layout(cfg);
%
%cd /analyse/4/Project0030/figures
%cfg        = [];
%cfg.layout = lay;
%cfg.zparam = 'stat';
%cfg.xlim   = [8 12];
%cfg.zlim   = 'absmax';
%cfg.interplimits = 'headleft';
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastalpha');
%cfg.xlim   = [18 22];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastbeta');
%cfg.xlim   = [40 80];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastgamma');
%
%cfg.xlim         = 56+[-12 12];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastgammalow12');
%cfg.xlim         = [85 100];
%figure;topoplotER(cfg,stat);
%print(gcf,'-depsc2','lrpdoublecontrastgammahigh12');
%
