% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var')
  error('subjectname needs to be defined');
end
if ~exist('range', 'var'), range = 'low'; end
if ~exist('doplanar', 'var'), doplanar = true; end
subject = vismot_subjinfo(subjectname);
 
[freq] = vismot_spectral_prepost_long(subject,'range', range, 'doplanar', doplanar);
filename = fullfile(subject.pathname,'freq',[subject.name,sprintf('freq_planar_long_%s', range)]);
save(filename, 'freq');
