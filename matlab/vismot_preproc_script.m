% this script performs preprocessing of a given subject, using pre-computed
% trl-matrices and pre-computed artifact definitions

if ~exist('subjectname', 'var')
  error('subjectname needs to be defined');
end
subject = vismot_subjinfo(subjectname);

for k = 1:5
  [~,trl(:,k)] = vismot_preproc_definetrial(subject,k);
end

alltrl = zeros(0,6);
for k = 1:size(trl,1)
  % if there is more than one run
  T = cat(1,trl{k,:});
  [srt,srtidx] = sort(T(:,1),'ascend');
  T = T(srtidx,:);
  alltrl = cat(1, alltrl, T);
end

% run it in a single shot, to get the gradiometer description consistent
% (in particular have the same denoising across conditions)
data = vismot_preproc_preprocessing(subject, alltrl);
  
for k = 1:5
  cfg = [];
  cfg.trials = find(data.trialinfo(:,1)==k);
  tmpdata    = ft_selectdata(cfg, data);
  
  switch k
    case 1
      data1 = tmpdata;
    case 2
      data2 = tmpdata;
    case 3
      data3 = tmpdata;
    case 4
      data4 = tmpdata;
    case 5
      data5 = tmpdata;
  end
  clear tmpdata;
end
filename = fullfile(subject.pathname,'data',[subject.name,'data']);
save(filename, 'data1', 'data2', 'data3', 'data4', 'data5');  
