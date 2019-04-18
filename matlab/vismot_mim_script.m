% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var')
	error('subjectname should be defined');
end
if ~exist('condition', 'var')
  condition = 1;
end
if numel(condition)==1
  split = true;
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
	[mim, parcellation] = vismot_mim_pre(subject,'prewhiten',true,'label',label,'split',split);
	filename = fullfile(subject.pathname,'mim',[subject.name,'_mim_pre']);
	
	mim = ft_struct2single(mim);
	save(filename, 'mim', 'parcellation');
	clear mim parcellation
end

