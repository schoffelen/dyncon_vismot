
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

subject = vismot_subjinfo(subjectname);

filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_enet.mat']);
load(filename, 'data', 'cfg')

cfg.design = cfg.design(randseq);

stat = ft_timelockstatistics(cfg, data);

filename = fullfile(subject.pathname,'pow', [subject.name, sprintf('_source3d4mm_pre_enet_rand%d.mat', randnr)]);
save(filename, 'stat')