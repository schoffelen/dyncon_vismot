datadir = '/project/3011085.03/analysis/source';
whichstat = 'statResp';%statCvsN statICvsN
% strings = {'statResp', 'statCvsN', 'statICvsN'};
% for l=1
%     whichstat = strings{l};
    frequency = [4:2:30 40:4:100];
    cnt = 0;
    for k = frequency
        fprintf('computing T-statistic for frequency %d Hz\n', k);
        
        d = dir(fullfile(datadir,sprintf('*3d4mm*pre_%03d.mat',k)));
        for m = 1:numel(d)
            dum = load(fullfile(d(m).folder,d(m).name),whichstat);
            tmp(m) = dum.(whichstat);
        end
        clear dumclear dum
        dat = cat(2,tmp.stat);
        n   = size(dat,2);
        for m = 1:n
            dat(:,m+n) = zeros(size((dat(:,m))));
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
source.gamma1 = nanmean(T(:,nearest(frequency, 56):nearest(frequency, 64)),2); % low gamma 50-70 Hz (56-64 Hz with 8 Hz smoothing )
source.gamma2 = nanmean(T(:,nearest(frequency, 76):nearest(frequency, 84)),2);% high gamma 70-90 Hz (76-84 Hz with 8 Hz smoothing )
source.alpha  = T(:,nearest(frequency, nearest(frequency, 10))); % classical alpha band 
source.beta1  = T(:,nearest(frequency, 16)); % low beta 12-20 Hz (16 Hz with 4 Hz smoothing)
source.beta2  = nanmean(T(:,nearest(frequency, 24):nearest(frequency, 26)),2); %20-30 Hz (24-26 Hz with 4 Hz smoothing)
    
    cfgi = [];
    cfgi.parameter = {'alpha' 'beta1' 'beta2' 'gamma1' 'gamma2'};
    source_int = ft_sourceinterpolate(cfgi, source, mri);
    
    source_int.alpha(~isfinite(source_int.alpha))=0;
    source_int.beta1(~isfinite(source_int.beta))=0;
    source_int.beta2(~isfinite(source_int.beta))=0;
    source_int.gamma1(~isfinite(source_int.gamma1))=0;
    source_int.gamma2(~isfinite(source_int.gamma2))=0;
    
    cmap = flipud(brewermap(64,'RdBu'));
    cfgp=[];
    cfgp.maskparameter = 'mask';
    cfgp.funcolormap = cmap;
    cfgp.maskstyle  = 'colormix';
    % cfgp.method     = 'slice';
    cfgp.nslices    =  30;
    cfgp.slicerange =  [40 150];
    cfgp.opacitylim = [2 4];
    
    % cfgp.method = 'surface';
    
    cfgp.funcolorlim = [-4 4];
    cfgp.funparameter = 'alpha';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('alpha');
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
%     h=(gcf); savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
    cfgp.funcolorlim = [-4 4];
    cfgp.funparameter = 'beta1';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('beta');
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
%     h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
        cfgp.funcolorlim = [-4 4];
    cfgp.funparameter = 'beta2';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('beta');
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
%     h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
    cfgp.funcolorlim = [-2.5 2.5];
    cfgp.method     = 'slice';
    cfgp.funparameter  = 'gamma1';
    ft_sourceplot(cfgp, source_int);
    title('gamma1');
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
%     h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
    cfgp.funcolorlim = [-2.5 2.5];
    cfgp.funparameter  = 'gamma2';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('gamma2'); 
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
%     h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
% end