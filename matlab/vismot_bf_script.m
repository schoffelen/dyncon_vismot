% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('frequency', 'var'),
  error('frequency should be defined');
end

if ~exist('subjectname', 'var'),
  error('subjectname needs to be defined');
end
subject = vismot_subjinfo(subjectname);
 
%[source, stat13, stat42] = vismot_bf_post(subject,'frequency',frequency);
%filename = fullfile(subject.pathname,'source',[subject.name,'source_post_',num2str(frequency)]);

load(fullfile(subject.pathname,'mri',sprintf('%s_sourcemodel3d6mm',subject.name)),'sourcemodel');
[source, stat13, stat42] = vismot_bf_post(subject,'frequency',frequency,'sourcemodel',sourcemodel);
filename = fullfile(subject.pathname,'source',[subject.name,'source3d_post_',num2str(frequency)]);

% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
  try,source(k).cfg = rmfield(source(k).cfg, {'grid' 'headmodel'}); end
  try,source(k).cfg.callinfo.usercfg = rmfield(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'}); end
end

%save(filename, 'source', 'stat13', 'stat42');  
save(filename, 'stat13', 'stat42');  

clear source stat13 stat42
