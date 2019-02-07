function visualize_cluster(clus, dim, inside, varargin)

%VISUALIZE_CLUSTER is a visualization tool that plots glass-brain projections
%of the clusters specified in the input
%
%Use as visualize_cluster(clus,dim,inside,hl,val,div,flag)
%
%clus can be a cell-array of structure-arrays, or a structure-array
%if it is a cell-array, the first structure of each of the cells is taken
%
%dim    (required) is the dimensionality of the volumes
%inside (required) the index vector to the inside voxels
%hl     (optional) is a matrix or cell-array containing the voxel indices
%val    (optional) is a matrix or cell-array containing information about the
%                     clusters (such as threshold peak-value etc FIXME
%div    (optional) is a scalar, a 1x2 vector, or an array specifying the number
%                     of clusters, the subdivision of the figure panel, or the
%                     list of clusters to be plotted

hl    = keyval('highlight', varargin);
val   = keyval('value',     varargin);
dat   = keyval('data',      varargin);
div   = keyval('div',       varargin);
flag  = keyval('flag',      varargin);
plotval = keyval('plotval', varargin);
wbg   = keyval('whitebg',   varargin);

if isempty(flag), flag = 0;       end
if isempty(plotval), plotval = 1; end
if isempty(wbg), wbg = 0; end

if iscell(clus) && isstr(clus{1}),
  fname    = clus{1};
  n1       = clus{2};
  n2       = clus{3};
  tmpclus1 = struct('vox',[],'val',[],'hit',[],'nvox',[],'thr',[]);
  tmpclus2 = struct('vox',[],'val',[],'hit',[],'nvox',[],'thr',[]);
  tmpclus1 = repmat(tmpclus1,[1 20]);
  tmpclus2 = repmat(tmpclus2,[1 20]);
  hl = [];
  %hack with some hard coded stuff
  load /home/jan/projects/ccc/matlab/grid08;
  for k = 1:20
    fprintf('loading in data\n');
    load([fname,'_',num2str(k)]);
    tmpclus1(k) = clus{n1}(1);
    tmpclus2(k) = cluss{n2}(1);
    %FIXME this has to be changed, do analyze_cluster2 convention
    ix                        = [];
    ia                        = find(ismember(grid.pos,location.pos,'rows'));
    [ix(:,1),ix(:,2),ix(:,3)] = ind2sub(dim,ia);
    hl{k}                     = ix;
    [d, state, num, tmpval]   = analyze_cluster(tmpclus1(k), dim, inside, ix, 1);
    val{k}                    = [tmpval; tmpclus1(k).thr; d];
    [d, state, num, tmpval]   = analyze_cluster(tmpclus2(k), dim, inside, ix, 1);
    vals{k}                   = [tmpval; tmpclus2(k).thr; d];
  end
  visualize_cluster(tmpclus1, dim, inside, hl, val,  div, flag);
  visualize_cluster(tmpclus2, dim, inside, hl, vals, div, flag);
elseif iscell(clus),
  if ~isfield(clus{1}(1),'state')
    tmpclus = struct('vox',[],'val',[],'hit',[],'nvox',[],'thr',[]);
  else
    tmpclus = struct('vox',[],'val',[],'hit',[],'nvox',[],'thr',[],...
    'd',[],'state',[],'num',[],'mindist',[]);
  end
  tmpclus = repmat(tmpclus,[1 length(clus)]);
  for k = 1:length(clus)
    if ~isempty(clus{k}),
      if isempty(hl),
        tmpclus(k) = clus{k}(1);
      else
        [d, state, num, tmpval] = analyze_cluster2(clus{k}, dim, inside, hl, 1); %FIXME the 1 is not OK always
        seld = find(d==min(d));
        tmpclus(k) = clus{k}(seld(1));
      end
    end
  end
  visualize_cluster(tmpclus, dim, inside, 'highlight', hl, 'value', val, 'div', div, 'flag', flag);
else
  if ~isempty(div),
    if prod(size(div))==1,
      selclus = 1:div; div = [];
    elseif length(div)==2,
      selclus = 1:min(length(clus), prod(div));
    elseif length(div)>2,
      selclus = div; div = [];
    end
  else
    selclus = 1:length(clus);
  end
  nclus = length(selclus);
  if isempty(div), div = [ceil(sqrt(nclus)) round(sqrt(nclus))]; end
  
  head = zeros(dim);
  head(inside) = 1;
  head1=-1.*double(squeeze(sum(head,1))==0)';
  head2=-1.*double(squeeze(sum(head,2))==0)';
  head3=-1.*double(squeeze(sum(head,3))==0)';
  
  siz1 = size(head1);
  siz2 = size(head2);
  siz3 = size(head3);
  head = zeros(siz2(1)+siz3(1)+4, siz2(2)+siz1(2)+4);
  head(2:siz2(1)+1,       2:siz2(2)+1)       = double(~head2);
  head(end-siz3(1):end-1, 2:siz3(2)+1)       = double(~head3);
  head(2:siz1(1)+1,       end-siz1(2):end-1) = double(~head1);
  
  blobs = zeros([size(head) nclus]);
  for k = 1:nclus
    ind = selclus(k);
    vol = zeros(dim);
    if ~isempty(clus(ind).vox)
      
      if plotval && isempty(dat),
        [ix(:,1),ix(:,2),ix(:,3),ix(:,4),ix(:,5),ix(:,6)] = ...
          ind2sub([dim dim], clus(ind).vox);
        for m = 1:size(ix,1)
          tmp = vol(ix(m,1),ix(m,2),ix(m,3));
          vol(ix(m,1),ix(m,2),ix(m,3)) = max(tmp, abs(clus(ind).val(m)));
          %        tmp = vol(ix(m,4),ix(m,5),ix(m,6));
          %        vol(ix(m,4),ix(m,5),ix(m,6)) = max(tmp, abs(clus(ind).val(m)));
        end
        sgn = sign(clus(ind).val(m));
        clear ix
      elseif plotval,
        
        if k==1, dummy=zeros(prod(dim)); dummy(inside,inside)=dat; end
        [ix(:,1),ix(:,2),ix(:,3),ix(:,4),ix(:,5),ix(:,6)] = ...
          ind2sub([dim dim], clus(ind).vox);
        for m = 1:size(ix,1)
          tmp = vol(ix(m,1),ix(m,2),ix(m,3));
          vol(ix(m,1),ix(m,2),ix(m,3)) = max(tmp, abs(dummy(clus(ind).vox(m))));
          tmp = vol(ix(m,4),ix(m,5),ix(m,6));
          vol(ix(m,4),ix(m,5),ix(m,6)) = max(tmp, abs(dummy(clus(ind).vox(m))));
        end
        sgn = sign(dummy(clus(ind).vox(m)));
        %keyboard
        clear ix;
        
      else
        [ix(:,1), ix(:,2)] = ind2sub([prod(dim) prod(dim)], clus(ind).vox);
        uix = unique(ix(:,2));
        for kk = 1:length(uix)
          dum = find(ix(:,1)==uix(kk));
          clus(ind).val(dum) = length(dum);
          %clus(ind).val(dum) = length(dum)./size(ix,1);
        end
        [ix(:,1),ix(:,2),ix(:,3),ix(:,4),ix(:,5),ix(:,6)] = ...
          ind2sub([dim dim], clus(ind).vox);
        for m = 1:size(ix,1)
          vol(ix(m,1),ix(m,2),ix(m,3)) = clus(ind).val(m);
        end
        %vol(ix(:,1)) = clus(ind).val;
        sgn = 1;
        clear ix
      end
    end
    
    dum = zeros(siz2(1)+siz3(1)+4, siz2(2)+siz1(2)+4);
    dum(2:siz2(1)+1,       2:siz2(2)+1)       = squeeze(max(vol,[],2))';
    dum(end-siz3(1):end-1, 2:siz3(2)+1)       = squeeze(max(vol,[],3))';
    dum(2:siz1(1)+1,       end-siz1(2):end-1) = squeeze(max(vol,[],1))';
    dum(dum==0) = nan;
    
    blobs(:,:,k) = sgn.*dum;
  end
  scale = max(abs(blobs(:)));
  scale = [-scale scale];
  siz   = size(blobs);
  
  ana   = zeros(div(1)*siz(1), div(2)*siz(2));
  fun   = zeros(div(1)*siz(1), div(2)*siz(2))*nan;
  for k = 1:nclus
    [ix,iy] = ind2sub(div, k);
    ana(end-ix*siz(1)+1:end-(ix-1)*siz(1),(iy-1)*siz(2)+1:iy*siz(2)) = head;
    fun(end-ix*siz(1)+1:end-(ix-1)*siz(1),(iy-1)*siz(2)+1:iy*siz(2)) = blobs(:,:,k);
  end
  
  ana = repmat(ana,[1 1 3]).*0.6;
  if wbg, ana(find(ana==0))=1; end;
  
  if flag==0,
    %plot in 1 figure window, which makes figure not possible to be saved with painters
    %as a renderer (to keep the vector objects), because rgb data and alpha masking is used
    figure;hold on;
    h  = imagesc(ana);
    hf = imagesc(fun);set(hf,'AlphaData',double(isfinite(fun)));
    caxis([scale]);
    axis xy; axis off; axis equal;
    for m=1:div(1)-1 h=plot([0.5 div(2)*siz(2)+0.5],[0.5 0.5]+m*siz(1)); set(h,'Color',ones(1,3)*(1-wbg));end
    for m=1:div(2)-1 h=plot([0.5 0.5]+m*siz(2),[0.5 div(1)*siz(1)+0.5]); set(h,'Color',ones(1,3)*(1-wbg));end
    
    if ~isempty(hl),
      if ~iscell(hl),
        hl = {hl};
        hl = repmat(hl, [1 nclus]);
      end
      for k = 1:nclus
        [ix,iy] = ind2sub(div, k);
        ix      = div(1)+1-ix;
        plot((iy-1)*siz(2)+hl{k}(:,1)+1,(ix-1)*siz(1)+hl{k}(:,3)+1,       'wo','markersize',4);
        plot((iy-1)*siz(2)+hl{k}(:,1)+1,(ix-1)*siz(1)+hl{k}(:,2)+3+dim(3),'wo','markersize',4);
        plot((iy-1)*siz(2)+hl{k}(:,2)+3+dim(1),(ix-1)*siz(1)+hl{k}(:,3)+1,'wo','markersize',4);
      end
    end
    
    if ~isempty(val),
      if ~iscell(val),
        val = {val};
        val = repmat(val, [1 nclus]);
      end
      for k = 1:nclus
        [ix,iy] = ind2sub(div, k);
        ix      = div(1)+1-ix;
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+6, num2str(val{k}(1),'%3.3f'),  'color', [1 1 1], 'fontsize', 6);
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+11, num2str(val{k}(2),'%3.4f'), 'color', [1 1 1], 'fontsize', 6);
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+16, num2str(val{k}(3),'%1.2f'), 'color', [1 1 1], 'fontsize', 6);
      end
    end
  else
    figure;hold on;
    h  = imagesc(ana);
    hf = imagesc(fun);set(hf,'AlphaData',double(isfinite(fun)));
    caxis([scale]);
    axis xy; axis off; axis equal;
    %set(gcf, 'PaperPositionMode', 'auto');
    set(gca, 'ActivePositionProperty', 'position');
    set(gca, 'Position', [0 0 1 1]);
    xlim1 = get(gca, 'xlim');
    ylim1 = get(gca, 'ylim');
    
    figure;hold on;
    set(gcf, 'InvertHardCopy', 'off');%, 'PaperPositionMode', 'auto');
    axis xy; axis off; axis equal;
    for m=1:div(1)-1 h=plot([0.5 div(2)*siz(2)+0.5],[0.5 0.5]+m*siz(1));set(h,'Color',ones(1,3)*(1-wbg)); end
    for m=1:div(2)-1 h=plot([0.5 0.5]+m*siz(2),[0.5 div(1)*siz(1)+0.5]);set(h,'Color',ones(1,3)*(1-wbg)); end
    
    if ~isempty(hl),
      if ~iscell(hl),
        hl = {hl};
        hl = repmat(hl, [1 nclus]);
      end
      for k = 1:nclus
        [ix,iy] = ind2sub(div, k);
        ix      = div(1)+1-ix;
        plot((iy-1)*siz(2)+hl{k}(:,1)+1,(ix-1)*siz(1)+hl{k}(:,3)+1,       'wo','markersize',4);
        plot((iy-1)*siz(2)+hl{k}(:,1)+1,(ix-1)*siz(1)+hl{k}(:,2)+3+dim(3),'wo','markersize',4);
        plot((iy-1)*siz(2)+hl{k}(:,2)+3+dim(1),(ix-1)*siz(1)+hl{k}(:,3)+1,'wo','markersize',4);
      end
    end
    
    if ~isempty(val),
      if ~iscell(val),
        val = {val};
        val = repmat(val, [1 nclus]);
      end
      for k = 1:nclus
        [ix,iy] = ind2sub(div, k);
        ix      = div(1)+1-ix;
      end
      for k = 1:nclus
        [ix,iy] = ind2sub(div, k);
        ix      = div(1)+1-ix;
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+6, num2str(val{k}(1),'%3.3f'),  'color', 'w', 'fontsize', 6);
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+11, num2str(val{k}(2),'%3.4f'), 'color', 'w', 'fontsize', 6);
        text((iy-1)*siz(2)+dim(1)+4, (ix-1)*siz(1)+dim(3)+16, num2str(val{k}(3),'%1.2f'), 'color', 'w', 'fontsize', 6);
      end
    end
    set(gca, 'ActivePositionProperty', 'position');
    set(gca, 'Position', [0 0 1 1]);
    set(gca, 'xlim', xlim1);
    set(gca, 'ylim', ylim1);
    xlim2 = get(gca, 'xlim');
    ylim2 = get(gca, 'ylim');
  end
end

