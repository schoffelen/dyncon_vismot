
load standard_sourcemodel3d4mm;
mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;
frequency = [10 22 38 42 58 62 78 82];
n=numel(list);

toi = 'pre';
roi_to = 'roi';
load(fullfile([alldir, 'analysis/source/roi.mat']));
nregions = 2*(size(roi,1)-1);

cnt=1;
for k=2:size(foi,1)
  tmpfreqs = foi{k,4};
  for subj=1:n
    sz = numel(tmpfreqs);
    tmpzx13 = [];
    tmpzx42 = [];
    tmpcoh1 = [];
    tmpcoh2 = [];
    tmpcoh3 = [];
    tmpcoh4 = [];
    for f = 1:sz
      
      filename = fullfile([datadir, sprintf('%s_coh6d4mm_%s_roi2%s_%03d', list{subj},toi,roi_to, tmpfreqs(f))]);
      dum = load(filename);
      
      tmpzx13(:,:,f) = dum.zx13;
      tmpzx42(:,:,f) = dum.zx42;
      tmpcoh1(:,:,f) = abs(dum.coh(1).coh);
      tmpcoh2(:,:,f) = abs(dum.coh(2).coh);
      tmpcoh3(:,:,f) = abs(dum.coh(3).coh);
      tmpcoh4(:,:,f) = abs(dum.coh(4).coh);
    end
    zx13{cnt}(subj,:,:) = nanmean(tmpzx13,3);
    zx42{cnt}(subj,:,:) = nanmean(tmpzx42,3);
    coh1{cnt}(subj,:,:) = nanmean(tmpcoh1,3);
    coh2{cnt}(subj,:,:) = nanmean(tmpcoh2,3);
    coh3{cnt}(subj,:,:) = nanmean(tmpcoh3,3);
    coh4{cnt}(subj,:,:) = nanmean(tmpcoh4,3);
  end
  cnt=cnt+1;
end

% % pool zx13 and zx42: flip zx42 and average
% for k=1:numel(zx42)
%   nloc = size(zx42{k},3)/2;
%   zx42

% calculate average within/between hemisphere coherence
% seperately for conditions 1-3 and 4-2
for k=1:numel(zx13)
  nloc = size(zx13{k},3)/2;
  x = ones(2*nloc);
  for subj=1:n
    x13tmp = squeeze(zx13{k}(subj,:,:));
    x13tmp(find(triu(x))) = nan;
    
    zx13_between(subj,k) = nanmean(nanmean(x13tmp(nloc+1:2*nloc, 1:nloc)));
    zx13_within(subj,k) = nanmean([nanmean(nanmean(x13tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x13tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    
    x42tmp = squeeze(zx42{k}(subj,:,:));
    x42tmp(find(triu(x))) = nan;
    zx42_between(subj,k) = nanmean(nanmean(x42tmp(nloc+1:2*nloc, 1:nloc)));
    zx42_within(subj,k) = nanmean([nanmean(nanmean(x42tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x42tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
    x1tmp = squeeze(coh1{k}(subj,:,:));
    x1tmp(find(triu(x))) = nan;
    coh1_between(subj,k) = nanmean(nanmean(x1tmp(nloc+1:2*nloc, 1:nloc)));
    coh1_within(subj,k) = nanmean([nanmean(nanmean(x1tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x1tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
        x2tmp = squeeze(coh2{k}(subj,:,:));
    x2tmp(find(triu(x))) = nan;
    coh2_between(subj,k) = nanmean(nanmean(x2tmp(nloc+1:2*nloc, 1:nloc)));
    coh2_within(subj,k) = nanmean([nanmean(nanmean(x2tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x2tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
        x3tmp = squeeze(coh3{k}(subj,:,:));
    x3tmp(find(triu(x))) = nan;
    coh3_between(subj,k) = nanmean(nanmean(x3tmp(nloc+1:2*nloc, 1:nloc)));
    coh3_within(subj,k) = nanmean([nanmean(nanmean(x3tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x3tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
        x4tmp = squeeze(coh4{k}(subj,:,:));
    x4tmp(find(triu(x))) = nan;
    coh4_between(subj,k) = nanmean(nanmean(x4tmp(nloc+1:2*nloc, 1:nloc)));
    coh4_within(subj,k) = nanmean([nanmean(nanmean(x4tmp(1:nloc, 1:nloc))),...
      nanmean(nanmean(x4tmp(nloc+1:2*nloc,nloc+1:2*nloc)))]);
    
  end
end

% pool conditions 13 and 42
zx_between = (zx13_between + zx42_between)./2;
zx_within = (zx13_within + zx42_within)./2;

% pool conditions 1+4 and 1+3
coh_between_C = (coh1_between + coh4_between)./2;
coh_between_IC = (coh3_between + coh2_between)./2;
coh_within_C = (coh1_within + coh4_within)./2;
coh_within_IC = (coh3_within + coh2_within)./2;

% average over frequencies
zx_between = nanmean(zx_between,2);
zx_within = nanmean(zx_within,2);
coh_between_C = nanmean(coh_between_C, 2);
coh_between_IC = nanmean(coh_between_IC, 2);
coh_within_C = nanmean(coh_within_C, 2);
coh_within_IC = nanmean(coh_within_IC, 2);


% we expect within hemisphere connectivity to be higher for conditions 1-3
% and 4-2, w.r.t. between hemispheres connectivity: i.e.
% zx_within>zxbetween
[w.H,w.P,w.CI,w.STATS] = ttest(zx_within, zx_between, 'tail', 'right');
within_vs_between_hemiC = w;

[w.H,w.P,w.CI,w.STATS] = ttest(coh_within_C, coh_within_IC, 'tail', 'right');
within_CvsIC = w;

[w.H,w.P,w.CI,w.STATS] = ttest(coh_between_C, coh_between_IC, 'tail', 'left');
between_CvsIC = w;


%% Is there any connectivity pattern significant individually (i.e. larger/smaller for C than IC)?
% combine 1-3 and 4-2
zx42_flipped = zx42;
for k=1:numel(zx42)
  nloc = size(zx13{k},3)/2;
  zx42_flipped{k}(:,1:nloc,1:nloc) = zx42{k}(:,nloc+1:2*nloc, nloc+1:2*nloc);
  zx42_flipped{k}(:,nloc+1:2*nloc, nloc+1:2*nloc) = zx42{k}(:,1:nloc,1:nloc);
  zx42_flipped{k}(:,nloc+1:2*nloc,1:nloc) = permute(zx42{k}(:,nloc+1:2*nloc,1:nloc), [1 3 2]);
  
  zx{k} = (zx13{k} + zx42_flipped{k})./2;
end
  
allconnectivity = cell(n,1);
for k=1:numel(zx)
  nloc = size(zx13{k},3)/2;
  x = ones(2*nloc);
  for subj=1:n
  tmpc = squeeze(zx{k}(subj,:,:));
  tmpc(find(triu(x))) = nan;
  tmpc = tmpc(:);
  tmpc = tmpc(~isnan(tmpc));
  allconnectivity{subj} = [allconnectivity{subj}, tmpc'];
  ncomparisson(k) = numel(tmpc);
  end
end
allconnectivity = cat(1, allconnectivity{:});

data = [];
data.label = {'coh'};
data.time = 1:size(allconnectivity,2);
data.trial = reshape(allconnectivity, [n 1 size(allconnectivity,2)]);
data.dimord = 'rpt_chan_time';

nul = data;
nul.trial = nul.trial*0;

cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.parameter = 'trial';
cfgs.alpha = 0.05;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.correctm = 'bonferroni';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';
stat = ft_timelockstatistics(cfgs, data, nul);

save(fullfile([alldir, 'analysis/stat_coh_pre.mat']), 'stat', 'allconnectivity', 'data', 'ncomparisson', 'zx', 'within_vs_between_hemiC', 'zx_within', 'zx_between', 'within_CvsIC', 'between_CvsIC','coh_within_C', 'coh_within_IC', 'coh_between_C', 'coh_between_IC');


%{
if ~exist('fliphemi', 'var'); fliphemi = false; end % not yet implemented. look at left and right hand resp trials seperately.
if ~exist('toi', 'var'); toi = 'pre'; end
if ~exist('include_neighb', 'var'); include_neighb = false; end
if ~exist('resamp', 'var'); resamp = false; end
if ~exist('spatsmooth', 'var'); spatsmooth=false; end
if ~exist('compute_var', 'var'); compute_var = 'Tstat'; end % can be 'avg'
if ~exist('doplot', 'var'); doplot=false; end
load standard_sourcemodel3d4mm

% prepare for hemiflipping right hand responses
insidepos = sourcemodel.pos(sourcemodel.inside==1,:);
[~, idx]=ismember(insidepos, sourcemodel.pos, 'rows');
s=sourcemodel;

mri = ft_read_mri('single_subj_T1_1mm.nii');
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

frequency = [10 22 38 42 58 62 78 82];
foi = [10 22 40 60 80];
n=19;
dat = zeros(74784, numel(frequency), n);
cnt = 1;
for k = frequency
    for m = 1:n
        filename = fullfile([datadir, sprintf('%s_coh6d4mm_%s_roi_%03d', list{m},toi, k)]);
        if include_neighb
            filename = fullfile([filename, '_neighb']);
        end
        if resamp
            filename = fullfile([filename, '_resamp']);
        end
        filename = fullfile([filename, '.mat']);
        dum = load(filename);
        zx13(m,cnt,:,:) = dum.zx13;
        if fliphemi
            % hemiflip right handed response trials
            s.x = zeros(size(s.inside,1), size(dum.zx42,2));
            s.x(idx,:) = dum.zx42;
            s=hemiflip(s, 'x');
            zx42(m,cnt,:,:) = s.x(s.inside==1,[2 1 4 3 6 5]);
        else
            zx42(m,cnt,:,:) = dum.zx42;
        end
    end
    clear dum
    cnt=cnt+1;
end

% average within frequency bands
[a1,a2,a3,a4] = size(zx42);
zx42tmp = zeros(a1,numel(foi), a3, a4);
zx42tmp(:,[1 2],:,:) = zx42(:,[1 2], :,:);
zx42tmp(:,3,:,:) = nanmean(zx42(:,[3 4],:,:),2);
zx42tmp(:,4,:,:) = nanmean(zx42(:,[5 6],:,:),2);
zx42tmp(:,5,:,:) = nanmean(zx42(:,[7 8],:,:),2);

zx13tmp = zeros(a1,numel(foi), a3, a4);
zx13tmp(:,[1 2],:,:) = zx13(:,[1 2], :,:);
zx13tmp(:,3,:,:) = nanmean(zx13(:,[3 4],:,:),2);
zx13tmp(:,4,:,:) = nanmean(zx13(:,[5 6],:,:),2);
zx13tmp(:,5,:,:) = nanmean(zx13(:,[7 8],:,:),2);

% only consider roi interactions
tmp = 1:numel(sourcemodel.inside);
tmp(sourcemodel.inside~=1)=[];
load('/project/3011085.03/analysis/stat_bf_pre.mat', 'l', 'r', 'freq_idx')
[~, inside_idx_l] = ismember(l,tmp);
[~, inside_idx_r] = ismember(r,tmp);
zx42tmp = zx42tmp(:, :, [inside_idx_l; inside_idx_r],:);
zx13tmp = zx13tmp(:, :, [inside_idx_l; inside_idx_r],:);

% only roi interactions in specific fois

%

% average within frequency bands
zx13(:,2,:,:) = nanmean(zx13(:,2:3,:,:),2);
zx13(:,3,:,:) = zx13(:,4,:,:);
zx13(:,4,:,:) = nanmean(zx13(:,5:6,:,:),2);
zx13(:,5,:,:) = nanmean(zx13(:,7:8,:,:),2);
zx13(:,6:end,:,:)=[];

zx42(:,2,:,:) = nanmean(zx42(:,2:3,:,:),2);
zx42(:,3,:,:) = zx42(:,4,:,:);
zx42(:,4,:,:) = nanmean(zx42(:,5:6,:,:),2);
zx42(:,5,:,:) = nanmean(zx42(:,7:8,:,:),2);
zx42(:,6:end,:,:)=[];

%% get ROI indices
load('/project/3011085.03/analysis/stat_bf_pre.mat', 'l', 'r')


%% spatial smoothing: average over neighbours (before tstat)
if spatsmooth
    [allneighb, seed, resolution] = find_neighbors(insidepos, sourcemodel);
    allneighb = reshape(permute(allneighb,[3,1,2]), [size(allneighb,1)*size(allneighb,3),3]);
    neighb_refindx = nan(size(allneighb,1),1);
    
    for m = 1:size(allneighb,1)
        [~,neighb_refindx(m)] = min( sum((insidepos-allneighb(m,:)).^2,2) ); % find the index of each ROI in insidepos.
    end
    [neighb_refindx, n_neighbors] = revise_neighbors(neighb_refindx, insidepos, resolution);
    
    tmp1 = zx13;
    tmp2 = zeros(size(zx13));
    index=1;
    for m=1:numel(n_neighbors)
        tmp2(:,:,m,:) = nanmean(tmp1(:,:,neighb_refindx(index:index+n_neighbors(m)),:), 3);
        index = (index+n_neighbors(m))+1;
    end
    zx13 = tmp2;
    % manually set coherence at original refindx to zero (because of contrast).
    for m=1:numel(refindx)
        zx13(:,:,refindx(m),m) = 0;
    end
    
    tmp1 = zx42;
    tmp2 = zeros(size(zx42));
    index=1;
    for m=1:numel(n_neighbors)
        tmp2(:,:,m,:) = nanmean(tmp1(:,:,neighb_refindx(index:index+n_neighbors(m)),:), 3);
        index = (index+n_neighbors(m))+1;
    end
    zx42 = tmp2;
    % manually set coherence at original refindx to zero (because of contrast).
    for m=1:numel(refindx)
        zx42(:,:,refindx(m),m) = 0;
    end
end
%% Average
if strcmp(compute_var, 'avg')
    % average over subjects
    zx13_avg = squeeze(nanmean(zx13,1));
    zx42_avg = squeeze(nanmean(zx42,1));
    
    % take only ROIs
    zx13ref=zx13_avg;
    zx42ref=zx42_avg;
    zx13ref(:,setdiff(1:size(zx13ref,2), refindx),:)=[];
    zx42ref(:,setdiff(1:size(zx42ref,2), refindx),:)=[];
    
    % make first four items in 2nd/3rd dimension left, the others right hemisphere
    zx13ref=zx13ref(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
    zx42ref=zx42ref(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
    
    % average over A+neighb-B coherence and B+neighb-A coherence.
    for k=1:size(zx13ref,3)
        for m=1:size(zx13ref,3)
            zx13ref(:,k,m) = (zx13ref(:,k,m)+zx13ref(:,m,k))/2;
            zx13ref(:,m,k) = zx13ref(:,k,m);
            
            zx42ref(:,k,m) = (zx42ref(:,k,m)+zx42ref(:,m,k))/2;
            zx42ref(:,m,k) = zx42ref(:,k,m);
        end
    end
end


%% T-statistic
if strcmp(compute_var, 'Tstat')
    [a,b,c,d] = size(zx13);
    % make FT structures, one for every ROI
    for k=1:size(zx13,4)
        coh(k)=sourcemodel;
    end
    % fill FT structures
    for k=1:d
        coh(k).freq=foi;
        coh(k).dimord = 'rpt_pos_freq';
        coh(k).coh13 = zeros(a,c,numel(foi));
        coh(k).coh42 = zeros(a,c,numel(foi));
        coh(k).inside=ones(c,1);
        coh(k).pos = coh(k).pos(sourcemodel.inside==1,:);
        coh(k).coh13 = permute(squeeze(zx13(:,:,:,k)), [1 3 2]);
        coh(k).coh42 = permute(squeeze(zx42(:,:,:,k)), [1 3 2]);
    end
    
    nul=coh(1);
    nul.coh13(:)=0;
    nul.coh42(:)=0;
    
    
    % statistics
    cfgs=[];
    cfgs.method = 'montecarlo';
    cfgs.statistic = 'depsamplesT';
    cfgs.alpha = 0.05;
    cfgs.ivar = 1;
    cfgs.uvar = 2;
    cfgs.design = [ones(1,n) ones(1,n)*2;1:n 1:n];
    cfgs.correctm = 'no';
    cfgs.numrandomization = 1000;
    cfgs.clusteralpha = 0.05;
    cfgs.correcttail = 'prob';
    for k=1:d
        cfgs.parameter = 'coh13';
        stat13(k) = ft_sourcestatistics(cfgs, coh(k), nul);
        cfgs.parameter = 'coh42';
        stat42(k) = ft_sourcestatistics(cfgs, coh(k), nul);
    end
    
    for k=1:d
        s13(k,:,:) = stat13(k).stat;
        s13ref(k,:,:) = stat13(k).stat(refindx,:);
        s42(k,:,:) = stat42(k).stat;
        s42ref(k,:,:) = stat42(k).stat(refindx,:);
    end
    s13 = permute(s13, [3 2 1]);
    s42 = permute(s42, [3 2 1]);
   
    s13ref = permute(s13ref, [3 1 2]);
    s42ref = permute(s42ref, [3 1 2]);
    s13ref = s13ref(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
    s42ref = s42ref(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);

    % hemiswap right hand responses
    % s42 = s42(:,[2 1 4 3 6 5],[2 1 4 3 6 5]);
    % % average between conditions
    % s = (s13+s42)/2;
    % make first three items in 2nd/3rd dimension left, 5-6th items right hemisphere
    % zx=zx(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
    if include_neighb; wh = 'seedneighb'; else wh='seed'; end
    filename = sprintf('/project/3011085.03/analysis/source/coh_Tstat_%s', wh);
    if spatsmooth
        filename = fullfile([filename '_spatsmooth']);
    end
    if resamp
        filename = fullfile([filename, '_resamp']);
    end
    save(filename, 'stat13', 'stat42', 'refindx', 'roi', 'zx13', 'zx42')
end

if doplot
    if strcmp(compute_var, 'avg')
        coh13ref = zx13ref;
        coh42ref = zx42ref;
        coh13 = zx13_avg;
        coh42 = zx42_avg;
    elseif strcmp(compute_var, 'Tstat')
        coh13ref = s13ref;
        coh42ref = s42ref;
        coh13 = s13;
        coh42 = s42;
    end
    % plot coherence between ROIs

    for k=1:numel(foi)
        figure(k);
        tmp=[abs(squeeze(coh13ref(k,:,:))), abs(squeeze(coh42ref(k,:,:)));];
        tmp = max(tmp(:));
        
        subplot(1,2,1);
        imagesc(squeeze(coh13ref(k,:,:)));
        title('left hand response'); caxis([-tmp tmp]); colorbar; pause(0.01)
        subplot(1,2,2);
        imagesc(squeeze(coh42ref(k,:,:)));
        title('right hand response'); caxis([-tmp tmp]); colorbar; pause(0.01)
        
        suptitle(sprintf('%d Hz', foi(k))); pause(0.001);
        saveas(gcf, sprintf('coh_%s_roi_%d.png',compute_var, foi(k)));
    end
    
    % look at full topography from ROIs to full brain.
    for l=[1 2]
        if l==1
            cohtmp = coh13;
            resphand = 'left';
        elseif l==2
            cohtmp = coh42;
            resphand = 'right';
        end
        clear source
        for k=1:6
            for m=1:numel(foi)
                source(k,m)=sourcemodel;
            end
        end
        for k=1:6
            for m=1:numel(foi)
                source(k,m).coh = zeros(size(sourcemodel.inside));
                source(k,m).coh(idx,:) = cohtmp(m,:,k)';
                source(k,m).freq = foi(m);
            end
        end
        
        cfg=[];
        cfg.parameter = 'coh';
        for k=1:6
            for m=1:numel(foi)
                z_int(k,m) = ft_sourceinterpolate(cfg, source(k,m), mri);
            end
        end
        
        
        cmap = flipud(brewermap(64,'RdBu'));
        cfgp=[];
        cfgp.funcolormap = cmap;
        cfgp.maskstyle  = 'colormix';
        cfgp.method     = 'slice';
        cfgp.nslices    =  30;
        cfgp.slicerange =  [40 150];
        cfgp.funparameter = 'coh';
        cfgp.maskparameter = 'coh'; %mask
        cfgp.visible = 'off';
        region = {'occ', 'occ', 'par', 'par', 'mot', 'mot'};
        if strcmp(resphand, 'left')
            side = {'ipsi', 'contra', 'ipsi', 'contra', 'ipsi', 'contra'};
        elseif strcmp(resphand, 'right')
            side = {'contra', 'ipsi', 'contra', 'ipsi', 'contra', 'ipsi',};
        end
        for k=1:6
            for m=1:5
                ft_sourceplot(cfgp, z_int(k,m));        title(sprintf('%s_%s_%d.png',region{k}, side{k}, foi(m)), 'Interpreter', 'none')
                saveas(gcf, sprintf('coh_%s_%s_%s_%d_%sresp.png', compute_var, region{k}, side{k}, foi(m), resphand))
            end
        end
    end
end
%}