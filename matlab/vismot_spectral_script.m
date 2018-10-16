% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var'),
  error('subjectname needs to be defined');
end
subject = vismot_subjinfo(subjectname);
 
[freqpre, freqpst] = vismot_spectral_prepost(subject,'foilim',[0 60],'smoothing',4);
filename = fullfile(subject.pathname,'freq','mve',[subject.name,'freq_axial_prepst_0-60']);
% save(filename, 'freqpre', 'freqpst');  
% save(filename, 'freqpst');  

%{ 
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
%}

%% statistic
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.numrandomization = 0;
n1 = size(freqpst(1).powspctrm,1);
n2 = size(freqpst(3).powspctrm,1);
cfg.design = [ones(1,n1) 2*ones(1,n2)];
stat13 = ft_freqstatistics(cfg, freqpst(1), freqpst(3));

n1 = size(freqpst(4).powspctrm,1);
n2 = size(freqpst(2).powspctrm,1);
cfg.design = [ones(1,n1) 2*ones(1,n2)];
stat42 = ft_freqstatistics(cfg, freqpst(4), freqpst(2));
stat42 = lrflip(stat42);

save(filename, 'freqpst', 'stat13', 'stat42'); 