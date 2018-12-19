% Define frequency (bands) of interest by looking at the C-IC post cue
% contrast, pooled for left hand and right hand responses (after flipping
% right hand responses. 

datadir = '/project/3011085.03/analysis/source';
whichstat = 'statResp';

frequency = [4:2:30 32:4:100];
cnt = 0;
for k = frequency
    fprintf('computing T-statistic for frequency %d Hz\n', k);
    
    d = dir(fullfile(datadir,sprintf('*3d4mm*post_%03d.mat',k)));
    for m = 1:numel(d)
        dum = load(fullfile(d(m).folder,d(m).name),whichstat);
        tmp(m) = dum.(whichstat);
    end
    clear dumclear dum
    dat = cat(2,tmp.stat);
    n   = size(dat,2);
    for m = 1:n
        dat(:,m+n) = nanmean(dat(:,m));
    end
    
    design   = [ones(1,n) ones(1,n)*2;1:n 1:n];
    cfg.ivar = 1;
    cfg.uvar = 2;
    tmp      = ft_statfun_depsamplesT(cfg, dat, design);
    cnt      = cnt+1;
    T(:,cnt) = tmp.stat;
    clear tmp;
    
end

load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
source = sourcemodel;
%%
for freq = frequency

source.tmp = nanmean(T(:,nearest(frequency, freq(1)):nearest(frequency, freq(end))),2);
cfgi = [];
cfgi.parameter = {'tmp'};
source_int = ft_sourceinterpolate(cfgi, source, mri);

source_int.tmp(~isfinite(source_int.tmp))=0;

cmap = flipud(brewermap(64,'RdBu'));
cfgp=[];
cfgp.funcolormap = cmap;
cfgp.maskstyle  = 'colormix';
cfgp.method     = 'slice';
cfgp.nslices    =  30;
cfgp.slicerange =  [40 150];
cfgp.opacitylim = [2 4];
cfgp.funparameter = 'tmp';
cfgp.funcolorlim = [-4 4];
cfgp.maskparameter = 'mask'; %tmp
ft_sourceplot(cfgp, source_int);
title(sprintf('%d Hz', freq));
h=frame2im(getframe(gcf)); imwrite(h, sprintf('%dHz.png', freq), 'PNG'); pause(0.01);
end