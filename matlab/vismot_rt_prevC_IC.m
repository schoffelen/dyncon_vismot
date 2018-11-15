
%% Are subjects faster when the previous trial was congruent?
load list
for k=1:19
subjectname = list(k);
subject = vismot_subjinfo(subjectname);
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = rmfield(alldata, 'data5');
alldata = vismot_data_reorder(alldata, 'previous');
C{k} = cat(1,alldata.data1.trialinfo,alldata.data4.trialinfo);
IC{k} = cat(1,alldata.data2.trialinfo,alldata.data3.trialinfo);
%C{k} = [alldata.data1.trialinfo(:,3);alldata.data4.trialinfo(:,3)];
%IC{k} = [alldata.data2.trialinfo(:,3);alldata.data3.trialinfo(:,3)];
end

for k=1:19
    IC_avg(k) = mean(IC{k}(:,3));
    C_avg(k) = mean(C{k}(:,3));
end

[H,P,CI,STATS] = ttest(C_avg, IC_avg)

% I would say that the thing to test is whether there's an interaction for
% the RT contingent on the congruency of the current trial, conditioned on the previous
% trial:
for k = 1:19
  rt(k,1) = mean(C{k}(ismember(C{k}(:,1),[1 4]),3)); %previous C, current C
  rt(k,2) = mean(C{k}(ismember(C{k}(:,1),[2 3]),3)); %previous C, current IC
  rt(k,3) = mean(IC{k}(ismember(IC{k}(:,1),[1 4]),3)); %previous IC, current C
  rt(k,4) = mean(IC{k}(ismember(IC{k}(:,1),[2 3]),3)); %previous IC, current IC
end 

% just judging the averages: theres a behavioral benefit for the current
% trial to be congruent (columns 1 and 3)
mean(rt)



%% Are subjects faster when the previous trial was congruent vs incongruent, while both required same response?

% previous and current response left
load list
for k=1:19
subjectname = list(k);
subject = vismot_subjinfo(subjectname);
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = rmfield(alldata, {'data2', 'data4', 'data5'});
alldata = vismot_data_reorder(alldata, 'previous');
C{k} = [alldata.data1.trialinfo(:,3)];
IC{k} = [alldata.data3.trialinfo(:,3)];
end

for k=1:19
nc = numel(C{k});
nic = numel(IC{k});
n = min([nc nic]);
randc = randperm(nc);
randic = randperm(nic);
C{k} = C{k}(randc(1:n));
IC{k} = IC{k}(randic(1:n));
end

for k=1:19
    IC_avg_left(k) = mean(IC{k});
    C_avg_left(k) = mean(C{k});
end

% previous and current response right
for k=1:19
subjectname = list(k);
subject = vismot_subjinfo(subjectname);
alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
alldata = rmfield(alldata, {'data1', 'data3', 'data5'});
alldata = vismot_data_reorder(alldata, 'previous');
C{k} = [alldata.data2.trialinfo(:,3)];
IC{k} = [alldata.data4.trialinfo(:,3)];
end

for k=1:19
nc = numel(C{k});
nic = numel(IC{k});
n = min([nc nic]);
randc = randperm(nc);
randic = randperm(nic);
C{k} = C{k}(randc(1:n));
IC{k} = IC{k}(randic(1:n));
end

for k=1:19
    IC_avg_right(k) = mean(IC{k});
    C_avg_right(k) = mean(C{k});
end


% combine
C_avg = [C_avg_left, C_avg_right];
IC_avg = [IC_avg_left, IC_avg_right];
[H,P,CI,STATS] = ttest(C_avg, IC_avg)
