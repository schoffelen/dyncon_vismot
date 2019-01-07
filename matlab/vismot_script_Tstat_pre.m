load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

whichstat = 'statResp';
frequency = [4:2:30 34:4:94];%[6 10 16 24 26 58 62 78 82];%
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 0;
for k = frequency
  for m = 1:n
    d = fullfile([datadir, sprintf('%ssource3d4mm_pre_%03d.mat', list{m}, k)]);
    dum = load(d, whichstat);
    tmp(m) = dum.(whichstat);
  end
  dat(:, cnt+1,1:n) = cat(2,tmp.stat);
  clear tmp dum
  cnt=cnt+1;
end
foi = {'theta', 6, 6
        'alpha', 10, 10
        'beta1', 16, 16
        'beta2', [24 26], 25
        'gamma1', [58 62], 60
        'gamma2', [78 82], 80};
frequency = [6 10 16 25 60 80];
    
dat = permute(dat, [3,1,2]);
source = sourcemodel;
source.dimord = 'rpt_pos_freq';
source.freq = [6 10 16 25 60 80];
source.stat(:,:,1:3) = dat(:,:,1:3);
source.stat(:,:,4) = nanmean(dat(:,:,4:5),3);
source.stat(:,:,5) = nanmean(dat(:,:,6:7),3);
source.stat(:,:,6) = nanmean(dat(:,:,8:9),3);

source_avg = rmfield(source, {'stat'});
source.avg = squeeze(nanmean(source.stat,1));

% Don't look at theta or low beta.
source.stat(:,:,[1 3]) = [];
source.avg(:,[1 3]) = [];
source.freq([1 3]) = [];

% load ROIs
load('/project/3011085.03/analysis/source/roi.mat', 'ROI');
for k=1:3
    idx_left(k) = find_dipoleindex(source, ROI{k,2});
    idx_right(k) = find_dipoleindex(source, ROI{k,3});
end
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'stat';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'no';
cfgs.numrandomization = 1000;
cfgs.clusteralpha = 0.05;
cfgs.correcttail = 'prob';
source.inside(:)=0;
source.inside([idx_left, idx_right]) = 1;
nul=source;
nul.stat(:)=0;

stat = ft_sourcestatistics(cfgs, source, nul);

% plot 2nd level tstat against zero
source.avg = stat.stat;
figure;
% occipital
subplot(1,3,1);
y = [source.avg(idx_left(1),1), source.avg(idx_left(1),3); source.avg(idx_right(1),1), source.avg(idx_right(1),3)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1'}, 'location', 'northwest')
title('occipital')

% parietal
subplot(1,3,2);
y = [source.avg(idx_left(2),1), source.avg(idx_left(2),3) source.avg(idx_left(2),4); source.avg(idx_right(2),1), source.avg(idx_right(2),3) source.avg(idx_right(2),4)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'alpha', 'gamma1', 'gamma2'}, 'location', 'northwest')
title('parietal')

% motor
subplot(1,3,3);
y = [source.avg(idx_left(3),2), source.avg(idx_left(3),4); source.avg(idx_right(3),2),source.avg(idx_right(3),4)]; 
bar(y);
set(gca,'xticklabel',{'left', 'right'});
legend({'beta2', 'gamma2'}, 'location', 'northwest')
title('motor')


% contrast left with right hemispheric ROI
stat_left = source.stat(:,idx_left,:); % rep_roi(left)_foi
stat_right = source.stat(:,idx_right,:); % rep_roi(right)_foi
stat_left_min_right = stat_left-stat_right;

[H,P,CI,STATS] = ttest(stat_left_min_right);
H = squeeze(H);
P = squeeze(P);
CI = squeeze(CI);
STATS = squeeze(STATS);
sprintf('OCCIPITAL Alpha: H=%d, p=%d; Gamma1: H=%d, p=%d', H(1,1), P(1,1), H(1,3), P(1,3))
sprintf('PARIETAL Alpha: H=%d, p=%d; Gamma1: H=%d, p=%d; Gamma2: H=%d, p=%d', H(2,1), P(2,1), H(2,3), P(2,3), H(2,4), P(2,4))
sprintf('MOTOR Beta2: H=%d, p=%d; Gamma2: H=%d, p=%d', H(3,2), P(3,2), H(3,4), P(3,4))


figure;
% occipital
subplot(1,3,1);
y = [STATS.tstat(1,1,1), STATS.tstat(1,1,3)]; 
bar(y);
set(gca,'xticklabel',{'alpha', 'gamma1'});
title('occipital')

% parietal
subplot(1,3,2);
y = [STATS.tstat(1,2,1), STATS.tstat(1,2,3) STATS.tstat(1,2,4)];
bar(y);
set(gca,'xticklabel',{'alpha', 'gamma1', 'gamma2'});
title('parietal')

% motor
subplot(1,3,3);
y = [STATS.tstat(1,3,2), STATS.tstat(1,3,4)]; 
bar(y);
set(gca,'xticklabel',{'beta2', 'gamma2'});
title('motor')


%%
%{
datadir = '/project/3011085.03/analysis/source';
whichstat = 'statResp';%statCvsN statICvsN
% strings = {'statResp', 'statCvsN', 'statICvsN'};
% for l=1
%     whichstat = strings{l};
    frequency = [10 24 26 58 62 78 82]%[4:2:30 40:4:100];
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
% source.beta1  = T(:,nearest(frequency, 16)); % low beta 12-20 Hz (16 Hz with 4 Hz smoothing)
source.beta2  = nanmean(T(:,nearest(frequency, 24):nearest(frequency, 26)),2); %20-30 Hz (24-26 Hz with 4 Hz smoothing)
    
    cfgi = [];
    cfgi.parameter = {'alpha' 'beta2' 'gamma1' 'gamma2'};
    source_int = ft_sourceinterpolate(cfgi, source, mri);
    
    source_int.alpha(~isfinite(source_int.alpha))=0;
%     source_int.beta1(~isfinite(source_int.beta1))=0;
    source_int.beta2(~isfinite(source_int.beta2))=0;
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
    
%     cfgp.funcolorlim = [-4 4];
    cfgp.funparameter = 'alpha';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('alpha');
    h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
    h=(gcf); savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
%     cfgp.funcolorlim = [-4 4];
%     cfgp.funparameter = 'beta1';
%     cfgp.method     = 'slice';
%     ft_sourceplot(cfgp, source_int);
%     title('beta');
%     h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
%     cfgp.method = 'surface';
%     ft_sourceplot(cfgp, source_int);
%     h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
%         cfgp.funcolorlim = [-4 4];
    cfgp.funparameter = 'beta2';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('beta');
    h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
    h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
%     cfgp.funcolorlim = [-2.5 2.5];
    cfgp.method     = 'slice';
    cfgp.funparameter  = 'gamma1';
    ft_sourceplot(cfgp, source_int);
    title('gamma1');
    h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
    h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
    
%     cfgp.funcolorlim = [-2.5 2.5];
    cfgp.funparameter  = 'gamma2';
    cfgp.method     = 'slice';
    ft_sourceplot(cfgp, source_int);
    title('gamma2'); 
    h=frame2im(getframe(gcf)); imwrite(h, sprintf('%s_%s.png', whichstat, cfgp.funparameter), 'PNG'); pause(0.01);
    cfgp.method = 'surface';
    ft_sourceplot(cfgp, source_int);
    h=gcf; savefig(h, sprintf('%s_%s', whichstat, cfgp.funparameter)); pause(0.01);
% end
%}