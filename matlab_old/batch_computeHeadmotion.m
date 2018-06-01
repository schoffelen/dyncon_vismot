addpath /home/jan/matlab/warping
subjinfo;
load /home/jan/projects/visuomotor/matlab/hmtlist

for k = 17 
  subject = SUBJ(k);
  if 1,
    cd([subject.pathname,'data']);
    load([subject.name,'data']);
    if ~strcmp(data1.label(end),'rv'),
      data1 = computeHeadmotion(data1, list, subject);
      save([subject.pathname,'data/',subject.name,'datatmp.mat'],'data1');clear data1
      data2 = computeHeadmotion(data2, list, subject);
      save([subject.pathname,'data/',subject.name,'datatmp.mat'],'data2','-append');clear data2
      data3 = computeHeadmotion(data3, list, subject);
      save([subject.pathname,'data/',subject.name,'datatmp.mat'],'data3','-append');clear data3
      data4 = computeHeadmotion(data4, list, subject);
      save([subject.pathname,'data/',subject.name,'datatmp.mat'],'data4','-append');clear data4
      data5 = computeHeadmotion(data5, list, subject);
      save([subject.pathname,'data/',subject.name,'datatmp.mat'],'data5','-append');clear data5
    else
      %already done
    end
  end
end

