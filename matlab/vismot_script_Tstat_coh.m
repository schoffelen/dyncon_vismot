    load standard_sourcemodel3d4mm;
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
        d = fullfile([datadir, sprintf('%scoh6d4mm_roi_%03d_resamp.mat', list{m}, k)]);
        dum = load(d);
        zx13(m,cnt,:,:) = dum.zx13;
        zx42(m,cnt,:,:) = dum.zx42;
        % hemiflip right handed response trials
        s.x = zeros(size(s.inside,1), size(dum.zx42,2));
        s.x(idx,:) = dum.zx42;
        s=hemiflip(s, 'x');
        zx42(m,cnt,:,:) = s.x(s.inside==1,[2 1 4 3 6 5]);
      end
      clear dum
      cnt=cnt+1;
    end
% mean over left and right hand response trials
zx = (zx13+zx42)/2;

      % average within frequency bands
    zx(:,2,:,:) = nanmean(zx(:,2:3,:,:),2);
    zx(:,3,:,:) = zx(:,4,:,:);
    zx(:,4,:,:) = nanmean(zx(:,5:6,:,:),2);
    zx(:,5,:,:) = nanmean(zx(:,7:8,:,:),2);
    zx(:,6:end,:,:)=[];

    % get ROI indices
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
      refindx_all = nan(size(roi,1),1);
      for m = 1:size(roi,1)
        [~,refindx(m)] = min( sum((insidepos-roi(m,:)).^2,2) ); % find the index of each ROI in insidepos.
        [~,refindx_all(m)] = min( sum((sourcemodel.pos-roi(m,:)).^2,2) ); % find the index of each ROI in sourcemodel.pos.
      end

      % average over subjects
      zx_avg = squeeze(nanmean(zx,1));
      
% average over neighbours
allneighb = find_neighbors(insidepos, sourcemodel);
allneighb = reshape(permute(allneighb,[3,1,2]), [size(allneighb,1)*size(allneighb,3),3]);
neighb_refindx = nan(size(allneighb,1),1);
  
  for m = 1:size(allneighb,1)
    [~,neighb_refindx(m)] = min( sum((insidepos-allneighb(m,:)).^2,2) ); % find the index of each ROI in insidepos.
  end
  tmp = reshape(neighb_refindx, [7, size(insidepos,1)])';
  zx_neighb = zeros(size(zx_avg));
  for k=1:size(insidepos,1)
      zx_neighb(:,k,:) = nanmean(zx_avg(:,tmp(k,:),:),2);
  end
  for k=1:numel(refindx)
      zx_neighb(:,refindx(k),k)=0;
  end
  
  zx_avg=zx_neighb;
      
% take only ROIs
zxref=zx_avg;
zxref(:,setdiff(1:size(zxref,2), refindx),:)=[];

% make first four items in 2nd/3rd dimension left, the others right hemisphere
zxref=zxref(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);

% average over A+neighb-B coherence and B+neighb-A coherence.
for k=1:size(zxref,3)
    for m=1:size(zxref,3)
        zxref(:,k,m) = (zxref(:,k,m)+zxref(:,m,k))/2;
        zxref(:,m,k) = zxref(:,k,m);
    end
end
    

for k=1:numel(foi)
    figure(k);
    imagesc(squeeze(zxref(k,:,:)));
    title(sprintf('%d Hz', foi(k))); pause(0.001);
    colorbar;
    tmp=abs(squeeze(zxref(k,:,:)));
    tmp = max(tmp(:));
    caxis([-tmp tmp])
    saveas(gcf, sprintf('coh_avg_roi_%d.png',foi(m)));
end

% look at full topography.
for k=1:6
    for m=1:numel(foi)
        z(k,m)=sourcemodel;
    end
end
for k=1:6
    for m=1:numel(foi)
    z(k,m).coh = zeros(size(sourcemodel.inside));
    z(k,m).coh(idx,:) = zx_avg(m,:,k)';
%     z(k).dimord = 'freq_pos';
    z(k,m).freq = foi(m);
    end
end

cfg=[];
cfg.parameter = 'coh';
for k=1:6
    for m=1:numel(foi)
    z_int(k,m) = ft_sourceinterpolate(cfg, z(k,m), mri);
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
        saveas(gcf, sprintf('%s_%s_%d.png',region{mod(k+2,3)+1}, side{mod(k+1,2)+1}, foi(m)))
    end
end
%%
[a,~,b,~] = size(zx13);
% make FT structures, one for every ROI
for k=1:size(zx13,4)
    coh(k)=sourcemodel;
end
% fill FT structures
for k=1:size(zx13,4)
coh(k).freq=foi;
coh(k).dimord = 'rpt_pos_freq';
coh(k).coh13 = zeros(a,b,numel(foi));
coh(k).coh42 = zeros(a,b,numel(foi));
coh(k).inside=ones(size(zx13,3),1);
coh(k).pos = coh(k).pos(sourcemodel.inside==1,:);
% average within frequency bands
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
    
for k=1:size(zx13,4)
    s13(k,:,:) = stat13(k).stat(refindx,:);
    s42(k,:,:) = stat42(k).stat(refindx,:);
end
s13 = permute(s13, [3,1,2]);
s42 = permute(s42, [3,1,2]);

% hemiswap right hand responses
s42 = s42(:,[2 1 4 3 6 5],[2 1 4 3 6 5]);
% average between conditions
s = (s13+s42)/2;
% average over subjects
% make first three items in 2nd/3rd dimension left, 5-6th items right hemisphere
zx=zx(:,[find(roi(:,1)<0); find(roi(:,1)>0)],[find(roi(:,1)<0); find(roi(:,1)>0)]);
  
  
  