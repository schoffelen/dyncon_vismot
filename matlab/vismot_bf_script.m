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
if ~exist('prewhiten', 'var')
  prewhiten = false;
end
if ~exist('lambda', 'var')
  lambda = [];
end
if ~exist('nrand', 'var')
  nrand = 100;
end
if isempty(smoothing)
  if frequency < 30
    smoothing = 4;
  else
    smoothing = 8;
  end
end
subject = vismot_subjinfo(subjectname);

load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45, statCvsIC] = vismot_bf_post(subject,'frequency',frequency,'sourcemodel',sourcemodel,'prewhiten',prewhiten, 'lambda', lambda, 'nrand', nrand);
filename = fullfile(subject.pathname,'source',[subject.name,'source3d4mm_post_',num2str(frequency,'%03d')]);
if istrue(prewhiten)
  filename = [filename '_prewhitened'];
end

% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
  source(k).cfg = removefields(source(k).cfg, {'grid' 'headmodel'});
  source(k).cfg.callinfo.usercfg = removefields(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'});
end

% hemiflip right handed response/ right hemifield
if isfield(stat42, 'statsmooth')
  parameter = {'stat', 'statsmooth'};
else
  parameter = 'stat';
end
stat42 = hemiflip(stat42, parameter);
stat43 = hemiflip(stat43, parameter);
stat25 = hemiflip(stat25, parameter);
stat45 = hemiflip(stat45, parameter);

% treat as if everything is left handed response
statResp = stat13;
statResp.stat = (statResp.stat + stat42.stat)/2;

% treat as if cue presented in left hemifield
statHemi = stat12;
statHemi.stat = (statHemi.stat + stat43.stat)/2;

statCvsN       = stat15;
statCvsN.stat  = (statCvsN.stat + stat45.stat)/2;
statICvsN      = stat35;
statICvsN.stat = (statICvsN.stat + stat25.stat)/2;

%save(filename, 'source', 'stat13', 'stat42','stat12', 'stat43','stat15','stat25','stat35','stat45', 'statResp', 'statHemi', 'statCvsN', 'statICvsN', 'statCvsIC', 'smoothing');

stat = rmfield(stat13,'stat');
stat.stat13 = stat13.stat;
stat.stat42 = stat42.stat;
stat.stat12 = stat12.stat;
stat.stat43 = stat43.stat;
stat.stat15 = stat15.stat;
stat.stat25 = stat25.stat;
stat.stat35 = stat35.stat;
stat.stat45 = stat45.stat;
stat.statResp = statResp.stat;
stat.statHemi = statHemi.stat;
stat.statCvsN = statCvsN.stat;
stat.statICvsN = statICvsN.stat;
stat.statCvsIC = statCvsIC.stat;

save(filename, 'stat', 'smoothing', 'lambda');

