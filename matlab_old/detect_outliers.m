function [trl,label] = detect_outliers(subject, condition, run)

if nargin<3,
  run = 1;
end

%switch denoising on if necessary
if isfield(subject, 'denoise'),
  denoiseflag = 1;
end

%do semi-automatic artifact detection
cd(subject.pathname);
cd('artifact');
load([subject.name,'eog-run',subject.runnames{run}(1)]);
eog  = cfg.artfctdef.eog;
load([subject.name,'jump-run',subject.runnames{run}(1)]);
jump = cfg.artfctdef.jump;
cd('../trl');
load([subject.name,'trl-run',subject.runnames{run}(1),'cnd',num2str(condition)]);
cfg.artfctdef.eog  = eog;
cfg.artfctdef.jump = jump;

%read header
hdr = read_header(cfg.datafile);

%do a first round of reading in the data with a highpassfilter applied;
%this in order to select a clean subset of trials for the weight computations

%update trial specification and remove UACurrent and line noise 
cfg.channel         = {'MEG'};
cfg.padding         = 5;
cfg.dftfilter       = 'yes';
cfg.artfctdef.reject = 'partial';
cfg.artfctdef.minaccepttim = 0.25;
cfg.artfctdef.type  = {'eog' 'jump'};
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 1;

cfg                 = rejectartifact(cfg);
data                = preprocessing(cfg);
dataorig            = data;
cfg.channel         = {'MEGREF'};
refdata             = preprocessing(cfg);
cfg.dftfilter       = 'no';
cfg.channel         = {'UACurrent'};
UACdata             = preprocessing(cfg);
dataorig            = appenddata([], dataorig, UACdata);
refdata             = appenddata([], refdata , UACdata);
if denoiseflag,
  cfg                 = [];
  cfg.denoise.channel = channelselection({'MEG'}, dataorig.label); 
  cfg.denoise.refchannel = 'UACurrent';
  cfg.denoise.hilbert = 'yes';
  data    = preprocessing(cfg, dataorig);
  cfg.denoise.channel = channelselection({'MEGREF'}, refdata.label);
  refdata = preprocessing(cfg, refdata);
end

cfg = [];
cfg.method      = 'summary';
cfg.keepchannel = 'yes';
cfg.channel     = 'MEG';
fprintf('performing rejectvisual to exclude outliers in MEG-data\n');
data            = rejectvisual(cfg, data);
close all;
%cfg.channel     = {'MEGREFG'};
%fprintf('performing rejectvisual to exclude outliers in REF-data\n');
%refdata         = rejectvisual(cfg, refdata);
%close all;
%cfg.channel     = {'MEGREFA'};
%fprintf('performing rejectvisual to exclude outliers in REF-data\n');
%refdata         = rejectvisual(cfg, refdata);
%close all;
%cfg.channel     = {'MEGREF' '-MEGREFA' '-MEGREFG'};
%fprintf('performing rejectvisual to exclude outliers in REF-data\n');
%refdata         = rejectvisual(cfg, refdata);
trl   = findcfg(data.cfg, 'trl');
label = data.label;

