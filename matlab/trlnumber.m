function trlnumber(subject, conditions)
if ~exist('conditions', 'var'); conditions = 'current'; end
if ischar(subject)
	subject = vismot_subjinfo(subject);
end

alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = vismot_data_reorder(alldata, conditions);
alldata = ft_appenddata([], alldata.data1, alldata.data2, alldata.data3, alldata.data4, alldata.data5);

trialnumber = alldata.trialinfo(:,2);
trialinfo = alldata.trialinfo;

filename = sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subject.name, conditions);
save(filename, 'trialnumber', 'trialinfo')