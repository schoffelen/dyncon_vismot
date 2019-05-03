
load list
subjectname = list{2};
subject = vismot_subjinfo(subjectname);

load(sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subject.name, 'previous'));
ntrials = numel(trialnumber);

for k=1:ntrials
  leaveouttrial = trialnumber(k);
  
qsubfeval(@vismot_execute_pipeline, 'vismot_mim_script', subject.name, {'prewhiten',true},{'split',false}, {'conditions', [1:5]}, {'dobaseline', 0}, {'doL1out', 1}, {'leaveouttrial', leaveouttrial}, 'memreq', 8*1024^3, 'timreq', 900, 'batchid', sprintf('mim_%s_%d', subjectname, leaveouttrial));

end
qsubfeval(@vismot_execute_pipeline, 'vismot_mim_script', subject.name, {'prewhiten',true},{'split',false}, {'conditions', [1:5]}, {'dobaseline', 0}, {'doL1out', 0}, 'memreq', 10*1024^3, 'timreq', 800, 'batchid', sprintf('mim_%s_all', subjectname));
