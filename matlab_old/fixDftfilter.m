subjinfo;
for k = [1 2 3 8]
  subject = SUBJ(k);
  %if length(subject.runnames)==1,
    cd(subject.pathname);
    cd('data');
    load([subject.name,'data']);
    [data1,data2,data3,data4,data5] = denoise_reject(subject,data1,data2,data3,data4,data5);
    save([subject.name,'datanew'],'data1','data2','data3','data4','data5');
    clear data1 data2 data3 data4 data5
  %end
end
