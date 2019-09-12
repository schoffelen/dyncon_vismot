function trlnumber(subject, conditions)
if ~exist('conditions', 'var'); conditions = 'current'; end
if ischar(subject)
	subject = vismot_subjinfo(subject);
end

alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = vismot_data_reorder(alldata, conditions);
if strcmp(conditions, 'previous')
  freqpre      = vismot_spectral(subject,'output','pow','keeptrials','yes', 'conditions','previous','toi','pre', 'dobaseline', 0, 'doL1out', 0);
  alldata = ft_appendfreq([], freqpre(1), freqpre(2), freqpre(3), freqpre(4));
else
  alldata = ft_appenddata([], alldata.data1, alldata.data2, alldata.data3, alldata.data4, alldata.data5);
end

trialnumber = alldata.trialinfo(:,2);
trialinfo = alldata.trialinfo;

filename = sprintf('/project/3011085.03/analysis/trl/%s_trialnumber_%s', subject.name, conditions);
save(filename, 'trialnumber', 'trialinfo')