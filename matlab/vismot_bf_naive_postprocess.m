load list
numrandomization = 100;
cnt=1;
for k=[1 7];
  k
  subjectname = list{k};
  subject = vismot_subjinfo(subjectname);
  filename = fullfile(subject.pathname,'pow', [subject.name, '_source3d4mm_pre_enet.mat']);
  tmp{cnt} = load(filename, 'stat','model', 'stat_perfreq');
  
    for l=1:numrandomization
      l
      filename = fullfile(subject.pathname,'pow', [subject.name, sprintf('_source3d4mm_pre_enet_rand%d.mat', l)]);
      randtmp{cnt,l} = load(filename);
    end
  cnt = cnt+1;
end

for k=1:2
  for l=1:numrandomization
    randacc(k,l) = randtmp{k,l}.stat.statistic.accuracy;
  end
  acc(k,1) = tmp{k}.stat.statistic.accuracy;
  model{k} = tmp{k}.model;
end




