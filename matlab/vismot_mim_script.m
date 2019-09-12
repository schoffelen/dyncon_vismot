% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var')
  error('subjectname should be defined');
end
if ~exist('condition', 'var')
  conditions = 1:5;
end
if ~exist('split', 'var')
  split = false;
end
if ~exist('dobaseline', 'var')
  dobaseline = 0;
end
if ~exist('doL1out', 'var')
  doL1out = false;
end
if ~exist('leaveouttrial', 'var')
  leaveouttrial = false;
end
subject = vismot_subjinfo(subjectname);

load(fullfile('/project/3011085.03/analysis/mri/Conte69_32k/atlas_subparc374_8k'));
label = atlas.parcellationlabel;
sel   = contains(label, '_1_');
sel   = sel|contains(label, '_2_');
sel   = sel|contains(label, '_3_');
sel   = sel|contains(label, '_4_');
sel   = sel|contains(label, '_5_');
sel   = sel|contains(label, '_6_');
sel   = sel|contains(label, '_7_');
sel   = sel|contains(label, '_17_');
sel   = sel|contains(label, '_18_');
sel   = sel|contains(label, '_19_');
label = label(sel);
if 1
  [mim, parcellation] = vismot_mim_pre(subject,'prewhiten',true,'label',label,'split',split, 'conditions', conditions, 'dobaseline', dobaseline, 'doL1out', doL1out, 'leaveouttrial', leaveouttrial);
  mim = ft_struct2single(mim);
  
  datadir = '/project_ext/3010029/reproducescript/analysis/';
  if doL1out
        filename = fullfile(datadir,'mim','tmp','singletrial', subject.name, [subject.name,'_mim_pre']);
%     filename = fullfile(subject.pathname,'mim','singletrial', subject.name, [subject.name,'_mim_pre']);
  else
        filename = fullfile(datadir,'mim','tmp', subject.name, [subject.name,'_mim_pre']);
%     filename = fullfile(subject.pathname,'mim',[subject.name,'_mim_pre']);
  end
  if ~split
    filename = [filename, '_all'];
  end
  if dobaseline
    filename = [filename, '_baseline'];
  end
  if doL1out
    filename = [filename, sprintf('_%d', leaveouttrial)];
    parcellation = [];
  end
  save(filename, 'mim', 'parcellation');
  clear mim parcellation
end

