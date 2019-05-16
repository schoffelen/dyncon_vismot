for l=1
load list
subjectname = list{l};
subject = vismot_subjinfo(subjectname);

load(sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subject.name, 'previous'));
ntrials = numel(trialinfo(:,2));
% in case not all jobs completed successfully
%{
filename = ['/project_ext/3010029/reproducescript/analysis/mim/singletrial/', subjectname, '/', '*.mat'];
d = dir(filename);
for k = 1:numel(d)
k;
% m{k} = load([d(k).folder,'/', d(k).name]);
donetrials(k) = str2num(char(extractBetween(d(k).name,"all_", ".mat")));
end
numel(trialnumber)
trialnumber = trialnumber(~ismember(trialnumber, donetrials));
numel(trialnumber)

%}
for k=1:numel(trialnumber)
  leaveouttrial = trialnumber(k);
  
qsubfeval(@vismot_execute_pipeline, 'vismot_mim_script', subject.name, {'prewhiten',true},{'split',false}, {'conditions', [1:4]}, {'dobaseline', 0}, {'doL1out', 1}, {'leaveouttrial', leaveouttrial}, 'memreq', 8*1024^3, 'timreq', 3600, 'batchid', sprintf('mim_%s_%d', subjectname, leaveouttrial));

end
qsubfeval(@vismot_execute_pipeline, 'vismot_mim_script', subject.name, {'prewhiten',true},{'split',false}, {'conditions', [1:4]}, {'dobaseline', 0}, {'doL1out', 0}, 'memreq', 10*1024^3, 'timreq', 3600, 'batchid', sprintf('mim_%s_all', subjectname));

clear donetrials trialnumber
end