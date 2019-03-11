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

if 0
	[granger, parcellation] = vismot_granger_pre(subject,'prewhiten',true,'label',label,'split',split);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre']);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger', 'parcellation');
	clear granger parcellation
end

if 0
	% do time-reversed granger
	[granger, parcellation] = vismot_granger_pre(subject, 'prewhiten',true,'label',label,'reverseflag', true, 'split',split);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre_timereversed']);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger');
	clear granger parcellation
end

if 1
	% do time-reversed granger
	[granger, parcellation] = vismot_granger_pre(subject, 'prewhiten',true,'label',label, 'conditions', condition, 'split',split);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre_conditionprev',num2str(condition)]);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger');
	clear granger parcellation
end

