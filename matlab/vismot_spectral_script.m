% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var'),
  error('subjectname needs to be defined');
end
subject = vismot_subjinfo(subjectname);
 
[freqpre, freqpst] = vismot_spectral_prepost(subject,'foilim',[0 60],'smoothing',4);
filename = fullfile(subject.pathname,'freq',[subject.name,'freq_axial_prepst_0-60']);
save(filename, 'freqpre', 'freqpst');  
clear freqpre freqpst
 
[freqpre, freqpst] = vismot_spectral_prepost(subject,'foilim',[40 120],'smoothing',8);
filename = fullfile(subject.pathname,'freq',[subject.name,'freq_axial_prepst_40-100']);
save(filename, 'freqpre', 'freqpst');  

[freqpre, freqpst] = vismot_spectral_prepost(subject,'foilim',[0 60],'smoothing',4,'doplanar',1);
filename = fullfile(subject.pathname,'freq',[subject.name,'freq_planar_prepst_0-60']);
save(filename, 'freqpre', 'freqpst');  
clear freqpre freqpst

[freqpre, freqpst] = vismot_spectral_prepost(subject,'foilim',[40 120],'smoothing',8,'doplanar',1);
filename = fullfile(subject.pathname,'freq',[subject.name,'freq_planar_prepst_40-100']);
save(filename, 'freqpre', 'freqpst');  
