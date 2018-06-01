% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition

if ~exist('subjectname', 'var'),
	error('subjectname should be defined');
end
if ~exist('condition', 'var'),
  condition = 1;
end
subject = vismot_subjinfo(subjectname);

if 0,
	[granger, parcellation] = vismot_granger_pre(subject);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre']);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger', 'parcellation');
	clear granger parcellation
end
if 0,
	% do time-reversed granger
	[granger, parcellation] = vismot_granger_pre(subject, 1);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre_timereversed']);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger', 'parcellation');
	clear granger parcellation
end
if 1,
	[granger, parcellation] = vismot_granger_pre(subject,0,1,condition);
	filename = fullfile(subject.pathname,'granger',[subject.name,'_granger_pre_conditionprev',num2str(condition)]);
	
	granger = ft_struct2single(granger);
	save(filename, 'granger');
	clear granger parcellation
end


