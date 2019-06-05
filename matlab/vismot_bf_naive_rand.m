
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/')
addpath('/project/3011085.03/scripts/fieldtrip/external/dmlt/external/svm/')

subject = vismot_subjinfo(subjectname);

filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_naive.mat']);
load(filename, 'data', 'cfg')
cfg.design = cfg.design(randseq);

stat = ft_timelockstatistics(cfg, data);

% filename = fullfile(subject.pathname,'pow', 'rand', [subject.name, sprintf('_source3d4mm_pre_naive_rand%d.mat', randnr)]);
filename = sprintf('/project_ext/3010029/reproducescript/vismsot/analysis/pow/rand/%s_source3d4mm_pre_naive_rand%d.mat', subject.name, randnr);
save(filename, 'statistic')