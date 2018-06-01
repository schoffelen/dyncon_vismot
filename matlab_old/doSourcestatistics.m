%---doSourcestatistics
subjinfo
addpath /home/jan/matlab/spm2

fname = [homejanmatlab,'/mri/templategrid6mm'];        
load(fname, 'grid');
cd([SUBJ(1).pathname,'source']);
cnt = 0;
for k = 1:length(SUBJ)
    try,
        fname = [SUBJ(k).name,'stat10'];
        fprintf('%s\n',fname);
        load(fname);
	stat13.pos = grid.pos;
	stat42.pos = grid.pos;
	stat13.dimord = 'pos_freq';
	stat42.dimord = 'pos_freq';
	stat13.freq   = 10;
	stat42.freq   = 10;
	s13{1} = stat13;
	s42{1} = stat42;
        fname = [SUBJ(k).name,'stat20'];
        fprintf('%s\n',fname);
        load(fname);
	stat13.pos = grid.pos;
	stat42.pos = grid.pos;
	stat13.dimord = 'pos_freq';
	stat42.dimord = 'pos_freq';
	stat13.freq   = 20;
	stat42.freq   = 20;
	s13{2} = stat13;
	s42{2} = stat42;
        fname = [SUBJ(k).name,'stat56'];
        fprintf('%s\n',fname);
        load(fname);
	stat13.pos = grid.pos;
	stat42.pos = grid.pos;
	stat13.dimord = 'pos_freq';
	stat42.dimord = 'pos_freq';
	stat13.freq   = 56;
	stat42.freq   = 56;
	s13{3} = stat13;
	s42{3} = stat42;

        s13 = selectdata(s13{:}, 'param', {'stat' 'stat2'});
        s42 = selectdata(s42{:}, 'param', {'stat' 'stat2'});
	cnt = cnt+1;
        allstat13{1,cnt} = s13;
        allstat42{1,cnt} = s42;
	clear s13 s42;
    end
end

gavg13 = selectdata(allstat13{:}, 'param', 'stat2');
gavg42 = selectdata(allstat42{:}, 'param', 'stat2');
Nsubj  = gavg13.dim(1);
clear allstat13 allstat42;

gavg13.stat2(Nsubj+1:2*Nsubj,:) = 0;
gavg42.stat2(Nsubj+1:2*Nsubj,:) = 0;
gavg13.dim(1) = 2*Nsubj;
gavg42.dim(1) = 2*Nsubj;

cfg           = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'pooledT';
cfg.parameter = 'stat2';
cfg.numrandomization = 500;
cfg.design    = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.ivar      = 1;
stat13        = sourcestatistics(cfg, gavg13);
stat13.dimord = 'pos_freq';
stat13.dim    = [size(stat13.pos,1) length(stat13.freq)];
stat42        = sourcestatistics(cfg, gavg42);
stat42.dimord = 'pos_freq';
stat42.dim    = [size(stat42.pos,1) length(stat42.freq)];

%---with the flipdim
%flipdim left-effector
dim  = pos2dim3d(grid.pos);
dum  = zeros(dim);
gavg = gavg42;
gavg.stat2(:,gavg.outside,:)   = 0;
gavg13.stat2(:,gavg.outside,:) = 0;
for k = 1:Nsubj
  for kk = 1:length(gavg13.freq)
    dum(:) = gavg13.stat2(k,:,kk);
    dum    = flipdim(dum,1);
    tmp    = reshape((gavg.stat2(k,:,kk)+dum(:)')./sqrt(2), dim);
    spm_smooth(tmp, tmp, 2);
    gavg.stat2(k,:,kk) = tmp(:);
  end
end
stat = sourcestatistics(cfg, gavg);
stat.dimord = 'pos_freq';
stat.dim    = [size(stat.pos,1) length(stat.freq)];

cfgp              = [];
cfgp.method       = 'ortho';
cfgp.funparameter = 'stat';
cfgp.interactive  = 'yes';
figure;sourceplot(cfgp, stat);


