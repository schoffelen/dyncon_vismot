% this script performs preprocessing of a given subject, using pre-computed
% trl-matrices and pre-computed artifact definitions

if ~exist('subjectname', 'var')
  error('subjectname needs to be defined');
end
subject = vismot_subjinfo(subjectname);

cfg = [];
cfg.dataset = fullfile(subject.rawpath,subject.name,subject.scanname,subject.sessionname,subject.emptyrun,subject.datafile);
cfg.channel = {'MEG' '-A40' '-A107' '-A157' '-A248'};
data = ft_preprocessing(cfg);

cfg = [];
cfg.length = 2;
cfg.overlap = 0.5;
data = ft_redefinetrial(cfg, data);

cfg = [];
cfg.demean = 'yes';
data = ft_preprocessing(cfg, data);

cfg = [];
cfg.method = 'summary';
data = ft_rejectvisual(cfg, data);
data.time(:) = data.time(1);

cfg          = [];
cfg.resamplefs = 300;
data         = ft_resampledata(cfg, data);

filename = fullfile(subject.pathname,'data',[subject.name,'emptyroom']);
save(filename, 'data');  
