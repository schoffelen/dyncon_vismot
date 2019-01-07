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
  
  load standard_sourcemodel3d4mm;
  insidepos = sourcemodel.pos(sourcemodel.inside,:);
  if islogical(sourcemodel.inside)
    insidevec = find(sourcemodel.inside);
  else
    insidevec = sourcemodel.inside;
  end
  refindx = nan(size(roi,1),1);
  for m = 1:size(roi,1)
    [~,refindx(m)] = min( sum((insidepos-roi(m,:)).^2,2) );
  end
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
if ~exist('smoothing', 'var')
    smoothing = 4;
end
subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[coh,zx13,zx42,looptime] = vismot_bf_pre_coh_roi(subject,'sourcemodel',sourcemodel,'frequency',frequency,'smoothing',smoothing,'nrand',nrand, 'refindx', refindx);

filename = fullfile(subject.pathname,'source',[subject.name,'coh6d4mm_roi_',num2str(frequency)]);
save(filename, 'zx13', 'zx42', 'coh');
