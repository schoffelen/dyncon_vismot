% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('frequency', 'var')
    error('frequency should be defined');
end

if ~exist('subjectname', 'var')
    error('subjectname needs to be defined');
end
if ~exist('smoothing', 'var')
    smoothing = 4;
end
subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d8mm',subject.name)),'sourcemodel');
[coh,zx13,zx42,looptime] = vismot_bf_pre(subject,'sourcemodel',sourcemodel,'frequency',frequency,'smoothing',smoothing,'nrand',0);

filename = fullfile(subject.pathname,'source',[subject.name,'coh6d8mm_',num2str(frequency)]);
save(filename, 'zx13', 'zx42');
