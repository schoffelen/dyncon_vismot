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
if isempty(smoothing)
    if frequency < 30
        smoothing = 4;
    else
        smoothing = 8;
    end
end
subject = vismot_subjinfo(subjectname);

load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[source, stat13, stat42, stat12, stat43, stat15, stat25, stat35, stat45] = vismot_bf_post(subject,'frequency',frequency,'sourcemodel',sourcemodel);
filename = fullfile(subject.pathname,'source',[subject.name,'source3d4mm_post_',num2str(frequency,'%03d')]);

% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
    try,source(k).cfg = rmfield(source(k).cfg, {'grid' 'headmodel'}); end
    try,source(k).cfg.callinfo.usercfg = rmfield(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'}); end
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
% statResp.statsmooth = (statResp.statsmooth + stat42.statsmooth)/2;
% treat as if cue presented in left hemifield
statHemi = stat12;
statHemi.stat = (statHemi.stat + stat43.stat)/2;
% statHemi.statsmooth = (statHemi.statsmooth + stat43.statsmooth)/2;

statCvsN = stat15;
statCvsN.stat = (statCvsN.stat + stat45.stat)/2;
statICvsN = stat35;
statICvsN.stat = (statICvsN.stat + stat25.stat)/2;

save(filename, 'source', 'stat13', 'stat42','stat12', 'stat43','stat15','stat25','stat35','stat45', 'statResp', 'statHemi', 'statCvsN', 'statICvsN', 'smoothing');

