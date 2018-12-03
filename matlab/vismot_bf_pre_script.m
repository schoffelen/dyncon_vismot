% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('frequency', 'var')
  error('frequency should be defined');
end

if ~exist('subjectname', 'var')
  error('subjectname needs to be defined');
end
if ~exist('smoothing', 'var')
  smoothing = [];
end
if ~exist('conditions', 'var')
  conditions = 'previous';
end

subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[source, stat] = vismot_bf_pre(subject,'sourcemodel',sourcemodel,'frequency',frequency,'smoothing',smoothing, 'conditions', conditions);

filename = fullfile(subject.pathname,'source',[subject.name,'source3d4mm_pre_',num2str(frequency)]);

% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
  try,source(k).cfg = rmfield(source(k).cfg, {'grid' 'headmodel'}); end
  try,source(k).cfg.callinfo.usercfg = rmfield(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'}); end
end
stat13 = stat.stat13;
stat42 = stat.stat42;
stat12 = stat.stat12;
stat43 = stat.stat43;
stat15 = stat.stat15;
stat25 = stat.stat25;
stat35 = stat.stat35;
stat45 = stat.stat45;

% hemiflip right handed response/ right hemifield
if isfield(stat42, 'tri')
  if isfield(stat42, 'statsmooth')
    parameter = {'stat', 'statsmooth'};
  else
    parameter = 'stat';
  end
  stat42 = hemiflip(stat42, parameter);
  stat43 = hemiflip(stat43, parameter);
  stat25 = hemiflip(stat25, parameter);
  stat45 = hemiflip(stat45, parameter);
elseif isfield(stat42, 'dim')
  stat42.stat = reshape(flip(reshape(stat42.stat,stat42.dim),1),[],1);
  stat43.stat = reshape(flip(reshape(stat43.stat,stat43.dim),1),[],1);
  %stat42.statsmooth = reshape(flip(reshape(stat42.statsmooth,stat42.dim),1),[],1);
  %stat43.statsmooth = reshape(flip(reshape(stat42.statsmooth,stat42.dim),1),[],1);
end

% treat as if everything is left handed response / cue presented in left
% hemifield.
statResp = stat13;
statResp.stat = (statResp.stat + stat42.stat)/2;
%statResp.statsmooth = (statResp.statsmooth + stat42.statsmooth)/2;
% also for neutral conditions
statCvsN = stat15;
statCvsN.stat = (stat15.stat + stat45.stat)/2;
statICvsN = stat35;
statICvsN.stat = (stat35.stat + stat25.stat)/2;

statHemi = stat12;
statHemi.stat = (statHemi.stat + stat43.stat)/2;
%statHemi.statsmooth = (statHemi.statsmooth + stat43.statsmooth)/2;


save(filename, 'stat13', 'stat42','stat12', 'stat43', 'statResp', 'statHemi', 'statCvsN', 'statICvsN');
% save(filename, 'stat13', 'stat42');

clear source stat13 stat42

