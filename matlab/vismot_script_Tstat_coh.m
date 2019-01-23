if ~exist('fliphemi', 'var'); fliphemi = false; end % not yet implemented. look at left and right hand resp trials seperately.
if ~exist('toi', 'var'); toi = 'pre'; end
if ~exist('include_neighb', 'var'); include_neighb = true; end
if ~exist('resamp', 'var'); resamp = false; end
if ~exist('spatsmooth_preT', 'var'); spatsmooth_preT=false; end
if ~exist('spatsmooth_postT', 'var'); spatsmooth_postT=false; end
if ~exist('compute_var', 'var'); compute_var = 'T'; end % can be 'avg'
if ~exist('doplot', 'var'); doplot=false; end
load standard_sourcemodel3d4mm

% prepare for hemiflipping right hand responses
insidepos = sourcemodel.pos(sourcemodel.inside==1,:);
[~, idx]=ismember(insidepos, sourcemodel.pos, 'rows');
s=sourcemodel;

mri = ft_read_mri('single_subj_T1_1mm.nii');
vismot_subjinfo;
alldir = '/project/3011085.03/';
datadir = fullfile([alldir, 'analysis/source/']);
load list;

frequency = [10 24 26 42 58 62 78 82];
foi = [10 25 42 60 80];
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
load('/project/3011085.03/analysis/source/roi.mat');
if iscell(ROI)
    % explicitly remove the last ROI, frontal region.
    roi = zeros(size(ROI,1)*2,3);
    for m = 1:size(ROI,1)
        roi((m-1)*2+1,:) = ROI{m,2};
        roi((m-1)*2+2,:) = ROI{m,3};
    end
end
roi = roi./10; % assume that the values were in mm, convert to cm

load standard_sourcemodel3d4mm;
insidepos = sourcemodel.pos(sourcemodel.inside,:);
if islogical(sourcemodel.inside)
    insidevec = find(sourcemodel.inside);
else
    insidevec = sourcemodel.inside;
end
refindx = nan(size(roi,1),1);
for m = 1:size(roi,1)
    [~,refindx(m)] = min( sum((insidepos-roi(m,:)).^2,2) ); % find the index of each ROI in insidepos.
end

%% spatial smoothing: average over neighbours (before tstat)
if spatsmooth_preT
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
    zx42_avg = squeeze(nanmean(zx13,1));
    
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
if strcmp('compute_var', 'T')
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
        coh(k).coh13 = squeeze(coh(k).coh13);
        coh(k).coh42 = squeeze(coh(k).coh42);
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
    for k=1:size(zx13,4)
        cfgs.parameter = 'coh13';
        stat13(k) = ft_sourcestatistics(cfgs, coh(k), nul);
        cfgs.parameter = 'coh42';
        stat42(k) = ft_sourcestatistics(cfgs, coh(k), nul);
    end
    
    for k=1:d
        s13(k,:,:) = stat13(k).stat;
        s42(k,:,:) = stat42(k).stat;
    end
    s13 = permute(s13, [3 2 1]);
    s42 = permute(s42, [3 2 1]);
    
    %% spatial smoothing: average over neighbours (after tstat)
    if spatsmooth_postT
        [allneighb, seed, resolution] = find_neighbors(insidepos, sourcemodel);
        allneighb = reshape(permute(allneighb,[3,1,2]), [size(allneighb,1)*size(allneighb,3),3]);
        neighb_refindx = nan(size(allneighb,1),1);
        
        for m = 1:size(allneighb,1)
            [~,neighb_refindx(m)] = min( sum((insidepos-allneighb(m,:)).^2,2) ); % find the index of each ROI in insidepos.
        end
        [neighb_refindx, n_neighbors] = revise_neighbors(neighb_refindx, insidepos, resolution);
        
        tmp1 = s13;
        tmp2 = zeros(size(s13));
        index=1;
        for m=1:numel(n_neighbors)
            tmp2(:,m,:) = nanmean(tmp1(:,neighb_refindx(index:index+n_neighbors(m)),:), 2);
            index = (index+n_neighbors(m))+1;
        end
        s13 = tmp2;
        % manually set coherence at original refindx to zero (because of contrast).
        for m=1:numel(refindx)
            s13(:,refindx(m),m) = 0;
        end
        
        tmp1 = s42;
        tmp2 = zeros(size(s42));
        index=1;
        for m=1:numel(n_neighbors)
            tmp2(:,m,:) = nanmean(tmp1(:,neighb_refindx(index:index+n_neighbors(m)),:), 2);
            index = (index+n_neighbors(m))+1;
        end
        s42 = tmp2;
        % manually set coherence at original refindx to zero (because of contrast).
        for m=1:numel(refindx)
            s42(:,refindx(m),m) = 0;
        end
    end
    
    % hemiswap right hand responses
    % s42 = s42(:,[2 1 4 3 6 5],[2 1 4 3 6 5]);
    % % average between conditions
    % s = (s13+s42)/2;
    % make first three items in 2nd/3rd dimension left, 5-6th items right hemisphere
    % zx=zx(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
end

if doplot
    if strcmp(compute_var, 'avg')
        coh13ref = zx13ref;
        coh42ref = zx42ref;
        coh13 = zx13_avg;
        coh42 = zx42_avg;
    elseif strcmp(compute_var, 'T')
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
        title('left hand response'); colorbar; caxis([-tmp tmp])
        subplot(1,2,2);
        imagesc(squeeze(coh42ref(k,:,:)));
        title('right hand response'); colorbar; caxis([-tmp tmp])
        
        suptitle(sprintf('%d Hz', foi(k))); pause(0.001);
        saveas(gcf, sprintf('coh_%s_roi_%d.png',compute_var, foi(m)));
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
        region = {'occ', 'par', 'mot'};
        side = {'ipsi', 'contra'};
        for k=1:6
            for m=1:5
                ft_sourceplot(cfgp, z_int(k,m));        title(sprintf('%s_%s_%d.png',region{mod(k+2,3)+1}, side{mod(k+1,2)+1}, foi(m)), 'Interpreter', 'none')
                saveas(gcf, sprintf('coh_%s_%s_%s_%d_%sresp.png', compute_var, region{mod(k+2,3)+1}, side{mod(k+1,2)+1}, foi(m)), resphand)
            end
        end
    end
end