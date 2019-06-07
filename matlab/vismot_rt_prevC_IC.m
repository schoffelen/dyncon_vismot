%% load data
load list

for k=1:numel(list)
  k
  subjectname = list(k);
  subject = vismot_subjinfo(subjectname);
  data{k} = load(fullfile(subject.pathname,'data',[subject.name,'data']));
  dataprev{k} = vismot_data_reorder(data{k}, 'previous');
  data{k} = ft_appenddata([], data{k}.data1, data{k}.data2, data{k}.data3,data{k}.data4,data{k}.data5)
  data{k} = data{k}.trialinfo;
  dataprev{k} = ft_appenddata([], dataprev{k}.data1, dataprev{k}.data2, dataprev{k}.data3,dataprev{k}.data4,dataprev{k}.data5);
  dataprev{k} = dataprev{k}.trialinfo;
end

% keep only unique trials
for k=1:numel(list)
  data{k} = unique(data{k}, 'rows', 'stable');
  dataprev{k} = unique(dataprev{k}, 'rows', 'stable');
end

subjects = vismot_subjinfo;
% add RT to trialinfo
for k=1:numel(list)
  tmprt = subjects(k).rt';
  tmprt = tmprt(:);
  data{k}(:,end+1) = tmprt(data{k}(:,2));
  dataprev{k}(:, end+1) = tmprt(dataprev{k}(:,2));
end

%% Are subjects faster for congruent trials than for incongruent trials?
for k=1:numel(list)
  idxC = data{k}(:,1)==1 | data{k}(:,1)==4;
  idxIC = data{k}(:,1)==3 | data{k}(:,1)==2;
  avgC(k) = mean(data{k}(idxC,end));
  avgIC(k) = mean(data{k}(idxIC,end));
end

[H,P,CI,STATS] = ttest(avgC, avgIC)
perc = mean(avgC./avgIC-1)*100;


% split for left and right response trials
% left reponse
for k=1:numel(list)
  idxC = data{k}(:,1)==1;
  idxIC = data{k}(:,1)==3;
  avgC_l(k) = mean(data{k}(idxC,end));
  avgIC_l(k) = mean(data{k}(idxIC,end));
end
[H_l,P_l,CI_l,STATS_l] = ttest(avgC_l, avgIC_l)

% right response
for k=1:numel(list)
  idxC = data{k}(:,1)==4;
  idxIC = data{k}(:,1)==2;
  avgC_r(k) = mean(data{k}(idxC,end));
  avgIC_r(k) = mean(data{k}(idxIC,end));
end
[H_r,P_r,CI_r,STATS_r] = ttest(avgC_r, avgIC_r)

filename = [subjects(1).pathname, '/rt/', 'stat_simon.mat'];
save(filename, 'data', 'avgC', 'avgIC', 'avgC_r', 'avgIC_r', 'avgC_l', 'avgIC_l','H', 'P', 'CI', 'STATS','perc', 'H_l', 'P_l', 'CI_l', 'STATS_l', 'H_r', 'P_r', 'CI_r', 'STATS_r')
%% Are subjects faster when the current and previous trial are of the same condition?
for k=1:numel(list)
  % current congruent
  idxC_C = (dataprev{k}(:,1)==1 & dataprev{k}(:,4)==1) | (dataprev{k}(:,1)==4 & dataprev{k}(:,4)==4) |...
    (dataprev{k}(:,1)==1 & dataprev{k}(:,4)==4) | (dataprev{k}(:,1)==4 & dataprev{k}(:,4)==1);
  idxN_C = (dataprev{k}(:,1)==1 & dataprev{k}(:,4)==5) | (dataprev{k}(:,1)==4 & dataprev{k}(:,4)==5);
  idxIC_C = (dataprev{k}(:,1)==1 & dataprev{k}(:,4)==2) | (dataprev{k}(:,1)==1 & dataprev{k}(:,4)==3) |...
    (dataprev{k}(:,1)==4 & dataprev{k}(:,4)==2) | (dataprev{k}(:,1)==4 & dataprev{k}(:,4)==3);
  
  %current incongruent
  idxIC_IC = (dataprev{k}(:,1)==3 & dataprev{k}(:,4)==3) | (dataprev{k}(:,1)==2 & dataprev{k}(:,4)==2) |...
    (dataprev{k}(:,1)==3 & dataprev{k}(:,4)==2) | (dataprev{k}(:,1)==2 & dataprev{k}(:,4)==3);
  idxN_IC = (dataprev{k}(:,1)==3 & dataprev{k}(:,4)==5) | (dataprev{k}(:,1)==2 & dataprev{k}(:,4)==5);
  idxC_IC = (dataprev{k}(:,1)==3 & dataprev{k}(:,4)==1) | (dataprev{k}(:,1)==2 & dataprev{k}(:,4)==4) |...
    (dataprev{k}(:,1)==3 & dataprev{k}(:,4)==4) | (dataprev{k}(:,1)==2 & dataprev{k}(:,4)==1);
  
  avgC_C(k) = mean(dataprev{k}(idxC_C,end));
  avgN_C(k) = mean(dataprev{k}(idxN_C,end));
  avgIC_C(k) = mean(dataprev{k}(idxIC_C,end));
 
  avgIC_IC(k) = mean(dataprev{k}(idxIC_IC,end));
  avgN_IC(k) = mean(dataprev{k}(idxN_IC,end));
  avgC_IC(k) = mean(dataprev{k}(idxC_IC,end));
end
ncond1=2;
ncond2=3;
n = numel(list);

S=1:n;
S = repmat(S, [1 ncond1*ncond2]);
F1 = [1*ones(1,n*ncond2), 2*ones(1,n*ncond2)];
F2 = repmat([1*ones(1,n) 2*ones(1,n), 3*ones(1,n)], [1 2]);
Y = [avgC_C avgN_C avgIC_C avgC_IC avgN_IC avgIC_IC];
factnames = {'current', 'previous'};
stat = rm_anova2(Y, S, F1, F2, factnames)
p = stat{4,6};
F = stat{4,5};
df = stat{4,3};
  
filename = [subjects(1).pathname, '/rt/', 'stat_gratton.mat'];
save(filename, 'dataprev', 'avgC_C', 'avgC_IC', 'avgIC_C', 'avgIC_IC', 'avgN_C', 'avgN_IC', 'stat', 'p','F','df')
%{  
  %% Are subjects faster when the previous trial was congruent?
  load list
  for k=1:19
    subjectname = list(k);
    subject = vismot_subjinfo(subjectname);
    alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
    alldata = vismot_data_reorder(alldata, 'previous');
    C{k} = cat(1,alldata.data1.trialinfo,alldata.data4.trialinfo);
    IC{k} = cat(1,alldata.data2.trialinfo,alldata.data3.trialinfo);
    N{k} = alldata.data5.trialinfo;
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
    %   rt(k,5) = mean(IC{k}
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
  
  
  %% Are Incongruent trials faster when visual information has to cross from left to right, vs right to left hemispheres?
  % (right hemisphere to left hemisphere (condition 2) faster then the other
  % way around (condition 3)
  
  load list
  for k=1:19
    subjectname = list(k);
    subject = vismot_subjinfo(subjectname);
    alldata = load(fullfile(subject.pathname,'data',[subject.name,'data']));
    alldata = rmfield(alldata, 'data5');
    C{k} = cat(1,alldata.data1.trialinfo,alldata.data4.trialinfo);
    IC{k} = cat(1,alldata.data2.trialinfo,alldata.data3.trialinfo);
    %IC{k} = [alldata.data2.trialinfo(:,3);alldata.data3.trialinfo(:,3)];
  end
  for k=1:19
    con2(k) = mean(IC{k}(IC{k}(:,1)==2,3));
    con3(k) = mean(IC{k}(IC{k}(:,1)==3,3));
  end
  [H,P,CI,STATS] = ttest(con2,con3, 'tail', 'left');
  %}
  
  
