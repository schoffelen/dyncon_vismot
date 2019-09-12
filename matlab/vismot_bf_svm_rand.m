
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

subject = vismot_subjinfo(subjectname);

filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_svm1342.mat']);
load(filename)

cfg.design = cfg.design(randseq);

stat13 = ft_timelockstatistics(cfg, data31);
stat42 = ft_timelockstatistics(cfg, data24);

statistic13 = stat13.statistic;
statistic42 = stat42.statistic;

filename = sprintf('/project_ext/3010029/reproducescript/vismsot/analysis/pow/rand/%s_source3d4mm_pre_svm1342_rand%d.mat', subject.name, randnr);
save(filename, 'statistic13', 'statistic42')