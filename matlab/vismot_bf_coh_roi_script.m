if ~exist('toi', 'var'); toi='pre'; end
if ~exist('include_neighb', 'var'); include_neighb=false; end
if ~exist('roi', 'var') && ~exist('refindx', 'var')
  if ~exist('frequency', 'var')
    error('frequency should be defined');
  end
if ~exist('nrand', 'var')
  nrand = 0;
end
if ~exist('subjectname', 'var')
    error('subjectname needs to be defined');
end
if ~exist('conditions', 'var')
    error('conditions need to be specified (current, previous, or current_previous')
end
if ~exist('smoothing', 'var')
    smoothing = 4;
end
if ~exist('lambda', 'var')
    lambda = [];
end
if ~isempty('lambda') && lambda(end)~='%'
    error('lambda should be specified in percentages')
end
  
  % use the rois defined in the file, at present with hard coded path
  % (suboptimal)
  d = load('/project/3011085.03/analysis/source/roi.mat');
  if iscell(d.roi)
    for m = 1:size(d.roi,1)-1
      roi(m,:) = d.roi{m+1,3};
      roi(m+size(d.roi,1)-1,:) = d.roi{m+1,4};
    end
  end
  if strcmp(d.unit, 'mm')
    roi = roi./10; % assume that the values were in mm, convert to cm
  end
end
if ~exist('refindx', 'var')
  
  load standard_sourcemodel3d4mm;
  %add neighbors to roi.
  if include_neighb
    [roi, roi_orig, resolution] = find_neighbors(roi, sourcemodel);
    roi = reshape(permute(roi,[3,1,2]), [size(roi,1)*size(roi,3),3]);
  end
  
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
  
  if include_neighb
  % exclude neighbors multiple copies of same neighbors and non-direct
  % neighbors
  [refindx, n_neighbors, index_orig_seed] = revise_neighbors(refindx, insidepos, resolution);
  else
      n_neighbors = 0; % no neigbors
      index_orig_seed = 1:numel(refindx);
  end
  ref.refindx = refindx;
  ref.n_neighbors = n_neighbors;
  ref.index_orig_seed = index_orig_seed;
  % refindx = find_dipoleindex(sourcemodel, roi); % Does not work yet.
end


subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[coh,zx13,zx42,looptime] = vismot_bf_coh_roi(subject,'toi', toi,'conditions', conditions, 'sourcemodel',sourcemodel,'lambda', lambda, 'frequency',frequency,'smoothing',smoothing,'nrand',nrand, 'ref', ref, 'include_neighb', include_neighb);

filename = fullfile(subject.pathname,'source', [subject.name, '_coh6d4mm_', sprintf('%s_', toi), sprintf('roi2%s_', roi_to),sprintf('%03d', frequency)] );
if include_neighb
    filename = fullfile([filename, '_neighb']); 
end
if nrand>0
    filename = fullfile([filename, '_resamp']);
end
save(filename, 'zx13', 'zx42', 'coh');
