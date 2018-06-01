cd /analyse/4/Project0030/trl/;
subjinfo

cnt = 0;
for k = [1:3 5:20]
  cnt = cnt+1;
  subject = SUBJ(k);
  subject.name
  nrun    = length(subject.runnames);
  for m = 1:nrun
    load([subject.name,'trl-run',num2str(subject.runnames{m}(1)),'longtrials']);
    allcfg{1,m} = cfg;
  end
  [a,b,c,rs{1,cnt}] = analyzeRT(allcfg);
  clear cfg allcfg;
end

for k = 1:length(rs)
  trimm(k,:) = rs{k}.ntrim./1017.25;
end
d1 = trimm(:,1)-trimm(:,3);
d2 = trimm(:,4)-trimm(:,2);

for k = 1:length(rs)
  ddatlft(:,:,k) = rs{k}.ddatlft;
  ddatrgt(:,:,k) = rs{k}.ddatrgt;
end
