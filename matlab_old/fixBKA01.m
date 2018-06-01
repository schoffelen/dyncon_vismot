subjinfo;
subject = SUBJ(17);

cd(subject.rawpath);
cd(subject.name);
cd(subject.cohname);

cfg = [];
cfg.channel  = {'MEG' 'UACurrent'};
cfg.datafile = 'e,rfhp1.0Hz,COH';
cfg.detrend  = 'yes';
data         = preprocessing(cfg);

topo1 = data.trial{1}(1:248,:)*data.trial{1}(249,:)';
topo2 = data.trial{2}(1:248,:)*data.trial{2}(249,:)';
topo3 = data.trial{3}(1:248,:)*data.trial{3}(249,:)';
topo4 = data.trial{4}(1:248,:)*data.trial{4}(249,:)';
topo5 = data.trial{5}(1:248,:)*data.trial{5}(249,:)';

cfg        = [];
cfg.layout = '4D248.lay';
lay        = prepare_layout(cfg);
[a,b]      = match_str(lay.label, data.label);

figure;topoplot([], lay.pos(a,1),lay.pos(a,2),topo1(b),lay.label(a));
%coil1 = rpa;
figure;topoplot([], lay.pos(a,1),lay.pos(a,2),topo2(b),lay.label(a));
%coil1 = lpa;
figure;topoplot([], lay.pos(a,1),lay.pos(a,2),topo3(b),lay.label(a));
%coil1 = nas;
figure;topoplot([], lay.pos(a,1),lay.pos(a,2),topo4(b),lay.label(a));
%coil1 = inion;
figure;topoplot([], lay.pos(a,1),lay.pos(a,2),topo5(b),lay.label(a));
%coil1 = cz;



