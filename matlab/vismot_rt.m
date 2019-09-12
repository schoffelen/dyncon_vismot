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

filename = [subjects(1).pathname, 'stat_behavior_simon.mat'];
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
  
filename = [subjects(1).pathname, 'stat_behavior_gratton.mat'];
save(filename, 'dataprev', 'avgC_C', 'avgC_IC', 'avgIC_C', 'avgIC_IC', 'avgN_C', 'avgN_IC', 'stat', 'p','F','df')
