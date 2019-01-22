if ~exist('toi', 'var'); toi='pre'; end
if ~exist('include_neighb', 'var'); include_neighb=true; end
if ~exist('roi', 'var') && ~exist('refindx', 'var')
  % as Nx3 matrix, in the same units as the sourcemodel
  %error('roi or refindx should be defined')
  
  % use the rois defined in the file, at present with hard coded path
  % (suboptimal)
  load('/project/3011085.03/analysis/source/roi.mat');
  if iscell(ROI)
    roi = zeros(size(ROI,1)*2,3);
    for m = 1:size(ROI,1)
      roi((m-1)*2+1,:) = ROI{m,2};
      roi((m-1)*2+2,:) = ROI{m,3};
    end
  end
  roi = roi./10; % assume that the values were in mm, convert to cm
end
if ~exist('refindx', 'var')
  
  load standard_sourcemodel3d4mm;
  %add neighbors to roi.
  if include_neighb
    roi = find_neighbors(roi, sourcemodel);
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
  % FIXME Now nearest neighbors are used, even if they're not direct
  % neighbors.
% refindx = find_dipoleindex(sourcemodel, roi); % Does not work yet.
end

if ~exist('frequency', 'var')
    error('frequency should be defined');
end

if ~exist('nrand', 'var')
  nrand = 100;
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
subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[coh,zx13,zx42,looptime] = vismot_bf_coh_roi(subject,'toi', toi,'conditions', conditions, 'sourcemodel',sourcemodel,'frequency',frequency,'smoothing',smoothing,'nrand',nrand, 'refindx', refindx, 'include_neighb', include_neighb);
if nrand>0
    filename = fullfile(subject.pathname,'source', [subject.name, '_coh6d4mm_', sprintf('%s_', toi), 'roi_',sprintf('%03d', frequency), '_resamp'] );
else
    filename = fullfile(subject.pathname,'source', [subject.name, '_coh6d4mm_', sprintf('%s_', toi), 'roi_',sprintf('%03d', frequency)] );
end
save(filename, 'zx13', 'zx42', 'coh');
