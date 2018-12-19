if ~exist('roi', 'var') && ~exist('refindx', 'var')
  % as Nx3 matrix, in the same units as the sourcemodel
  error('roi or refindx should be defined')
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
