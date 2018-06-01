subjinfo;
cnt = 0;
for k = 1:length(SUBJ)
  %fname = [SUBJ(k).name,'lrp'];
  fname = [SUBJ(k).name,'tfr010'];
  try,
    load(fname, 'allfreq');
    cnt = cnt+1;
  end
  for j = 1:4
    dat{cnt,j}   = squeeze(allfreq(j).stat);
    label{cnt,j} = allfreq(j).label;
    time{cnt,j}  = allfreq(j).time;
  end
end

cfg.layout = '4D248.lay';
lay = prepare_layout(cfg);
ok  = zeros(248,1);
for k = 1:numel(label)
  [a,b] = match_str(lay.label, label{k});
  ok(a) = ok(a)+1;
end
oklabel = lay.label(ok==numel(label));

T = zeros(numel(dat),length(oklabel),length(time{1,1}));
n = size(dat,1);
for k = 1:size(dat,1)
  for m = 1:size(dat,2)
    [a,b] = match_str(oklabel, label{k,m});
    T((m-1)*n+k,:,:) = dat{k,m}(b,:);
    newlabel{k,m}    = label{k,m}(b);
  end
end

tlck        = [];
tlck.dimord = 'rpt_chan_time';
tlck.time   = time{1,1};
tlck.avg    = T;
tlck.label  = oklabel;

cfg                  = [];
cfg.method           = 'montecarlo';
cfg.numrandomization = 0;
cfg.statistic        = 'anova2x2rm';
cfg.f                = 'main1';
cfg.precondition     = 'after';
cfg.parameter        = 'avg';
cfg.design           = [repmat(1:n,[1 4]);
                        ones(1,n*2) ones(1,n*2)*2;...
                        ones(1,n)   ones(1,n)*2 ones(1,n) ones(1,n)*2];
cfg.uvar             = 1;
cfg.ivar             = [2 3];
stat                 = timelockstatistics(cfg, tlck);


cnt = zeros(248,1);
for k = 1:length(dat)
  [a,b] = match_str(lay.label, dat{k}.label);
  cnt(a) = cnt(a)+1;
end

sel     = find(cnt==length(dat));
channel = lay.label(sel);
for k = 1:length(dat)
  dat2{k} = selectdata(dat{k},'channel',channel);
end
dat    = dat2;clear dat2;
for k = 1:length(dat)
  dat{k}.powspctrm = dat{k}.avg;
  dat{k} = rmfield(dat{k},'avg');
  dat{k}.dimord = 'chan_freq';
  dat{k}.freq = dat{k}.time;
end  
alldat = selectdata(dat{:},'param','powspctrm');
nsubj  = size(alldat.powspctrm,1);
alldat.powspctrm(nsubj+1:2*nsubj,:,:) = 0;

cd /home/jan/projects/visuomotor/matlab
load leftchanlist;
cd /analyse/4/Project0030/data
load CDE04data data1;
alldat.grad          = data1.grad;
alldat               = rmfield(alldat,'time');
cfg                  = [];
cfg.grad             = data1.grad;
cfg.neighbourdist    = 0.037;
neighb               = neighbourselection(cfg);
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.numrandomization = 500;
cfg.statistic        = 'pooledT';
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.design = [ones(1,nsubj) ones(1,nsubj)*2; 1:nsubj 1:nsubj];
cfg.ivar   = 1;
cfg.uvar   = 2;
cfg.channel = channel;
cfg.neighbours = neighb;
stat = freqstatistics(cfg, alldat);

cfg = [];
cfg.layout = '4D248.lay';
lay = prepare_layout(cfg);

cd /analyse/4/Project0030/figures
cfg        = [];
cfg.layout = lay;
cfg.zparam = 'stat';
cfg.xlim   = [8 12];
cfg.zlim   = 'absmax';
cfg.interplimits = 'headleft';
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','lrpdoublecontrastalpha');
cfg.xlim   = [18 22];
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','lrpdoublecontrastbeta');
cfg.xlim   = [40 80];
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','lrpdoublecontrastgamma');

cfg.xlim         = 56+[-12 12];
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','lrpdoublecontrastgammalow12');
cfg.xlim         = [85 100];
figure;topoplotER(cfg,stat);
print(gcf,'-depsc2','lrpdoublecontrastgammahigh12');

