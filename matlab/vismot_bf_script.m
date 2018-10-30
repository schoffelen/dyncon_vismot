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

%[source, stat13, stat42, stat12, stat43] = vismot_bf_post(subject,'frequency',frequency, 'smoothing', smoothing);
%filename = fullfile(subject.pathname,'source','mve',[subject.name,'source_post_',num2str(frequency)]);
% filename = fullfile(subject.pathname,'source','mve','cmnfiltpost',[subject.name,'source_post_',num2str(frequency)]);
 

%load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d6mm',subject.name)),'sourcemodel');
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');

[source, stat13, stat42, stat12, stat43] = vismot_bf_post(subject,'frequency',frequency,'sourcemodel',sourcemodel);
%filename = fullfile(subject.pathname,'source',[subject.name,'source3d_post_',num2str(frequency)]);
filename = fullfile(subject.pathname,'source',[subject.name,'source3d4mm_post_',num2str(frequency)]);

% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
    try,source(k).cfg = rmfield(source(k).cfg, {'grid' 'headmodel'}); end
    try,source(k).cfg.callinfo.usercfg = rmfield(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'}); end
end

% hemiflip right handed response/ right hemifield
if isfield(stat42, 'tri')
  n = size(stat42.stat,1)./2;
  stat42.stat = stat42.stat([n+(1:n) 1:n],1);
  stat42.statsmooth = stat42.statsmooth([n+(1:n) 1:n],1);
  n = size(stat43.pos,1)./2;
  stat43.stat = stat43.stat([n+(1:n) 1:n],1);
  stat43.statsmooth = stat43.statsmooth([n+(1:n) 1:n],1);
elseif isfield(stat42, 'dim')
  stat42.stat = reshape(flip(reshape(stat42.stat,stat42.dim),1),[],1);
  stat43.stat = reshape(flip(reshape(stat42.stat,stat42.dim),1),[],1);
  stat42.statsmooth = reshape(flip(reshape(stat42.statsmooth,stat42.dim),1),[],1);
  stat43.statsmooth = reshape(flip(reshape(stat42.statsmooth,stat42.dim),1),[],1);
  
end

% treat as if everything is left handed response / cue presented in left
% hemifield.
statResp = stat13;
statResp.stat = (statResp.stat + stat42.stat)/2;
statResp.statsmooth = (statResp.statsmooth + stat42.statsmooth)/2;
statHemi = stat12;
statHemi.stat = (statHemi.stat + stat43.stat)/2;
statHemi.statsmooth = (statHemi.statsmooth + stat43.statsmooth)/2;

save(filename, 'stat13', 'stat42','stat12', 'stat43', 'statResp', 'statHemi');
% save(filename, 'stat13', 'stat42');

clear source stat13 stat42
