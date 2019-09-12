function [dataout] = vismot_data_reorder(datain, conditions)

switch conditions
  case 'current'
    % nothing needed
    dataout = datain;
  case 'previous'
    % order according to previous trial.
    datain=removefields(datain, 'trl');
    
    % sanity check on input
    if numel(fieldnames(datain))~=5
      warning('unexpected number of data arguments in input');
    end
    
    fd = fieldnames(datain);
    for k=1:length(fd)
        tmpdat{k} = datain.(fd{k});
    end
%     data = ft_appenddata([],datain.data1,datain.data2,datain.data3,datain.data4,datain.data5);
    data = ft_appenddata([],tmpdat{:});
    clear datain;
    
    % assign new condition number based on the condition of the previous
    % trial -> this is based on getting information from the data that is
    % artifact-free. It can be improved by using the original trl-info
    % FIXME
    c = data.trialinfo(:,1);
    t = data.trialinfo(:,2);
    
    newc = zeros(numel(data.trial),1);
    for k = 1:numel(c)
      sel = find(t==t(k)-1);
      if numel(sel)==1
        newc(k) = c(sel);
      end
    end
    % add the condition of the previous trial to the trialinfo.
    data.trialinfo(:,end+1) = newc;
    
    % remove trials that were not preceded by a different trial, but
    % instead were preceded by the baseline (trials were presented in
    % blocks of 4).
    cfg=[];
    trialnumber = data.trialinfo(:,2);
    not_first_trial_of_block = find(mod(trialnumber, 4)~=1);
    cfg.trials = not_first_trial_of_block;
    data = ft_selectdata(cfg, data);
    newc = newc(not_first_trial_of_block,:);
    
    cfg = [];
    cfg.trials    = find(newc==1);
    dataout.data1 = ft_selectdata(cfg, data);
    cfg.trials    = find(newc==2);
    dataout.data2 = ft_selectdata(cfg, data);
    cfg.trials    = find(newc==3);
    dataout.data3 = ft_selectdata(cfg, data);
    cfg.trials    = find(newc==4);
    dataout.data4 = ft_selectdata(cfg, data);
    cfg.trials    = find(newc==5);
    dataout.data5 = ft_selectdata(cfg, data);
    
  case 'current_previous'
    % order according to current trial, conditioned on previous
    trl = datain.trl';
    trl = trl(:); % this vector contains conditions of all trials in chronological order
    datain = removefields(datain, 'trl');
    % sanity check on input
    if numel(fieldnames(datain))~=5
      warning('unexpected number of data arguments in input');
    end
    
    fd = fieldnames(datain);
    for k=1:length(fd)
        tmpdat{k} = datain.(fd{k});
    end
%     data = ft_appenddata([],datain.data1,datain.data2,datain.data3,datain.data4,datain.data5);
    data = ft_appenddata([],tmpdat{:});
    clear datain;
    
    % assign new condition number based on the condition of the previous
    % trial -> this is based on getting information from the data that is
    % artifact-free. It can be improved by using the original trl-info
    % FIXME
    c = data.trialinfo(:,1);
    t = data.trialinfo(:,2);
    
    prevc = zeros(numel(data.trial),1);
    for k = 1:numel(c)
      if t(k)~=1
        prevc(k) = trl(t(k)-1);
      end
    end
    % add the condition of the previous trial to the trialinfo.
    data.trialinfo(:,end+1) = prevc;
    
    % remove trials that were not preceded by a different trial, but
    % instead were preceded by the baseline (trials were presented in
    % blocks of 4).
    cfg=[];
    trialnumber = data.trialinfo(:,2);
    not_first_trial_of_block = find(mod(trialnumber, 4)~=1);
    cfg.trials = not_first_trial_of_block;
    data = ft_selectdata(cfg, data);
    c = c(not_first_trial_of_block,:);
    prevc = prevc(not_first_trial_of_block,:);
    
    cfg=[];
    cfg.trials = find((c==1) & (prevc==1 | prevc==4)); % C-C
    dataout.data1L = ft_selectdata(cfg, data);
    cfg.trials = find((c==4) & (prevc==1 | prevc==4)); % C-C
    dataout.data1R = ft_selectdata(cfg, data);
    cfg.trials = find((c==1) & (prevc==2 | prevc==3)); % C-IC
    dataout.data2L = ft_selectdata(cfg, data);
    cfg.trials = find((c==4) & (prevc==2 | prevc==3)); % C-IC
    dataout.data2R = ft_selectdata(cfg, data);
    cfg.trials = find((c==1) & prevc==5); % C-N
    dataout.data3L = ft_selectdata(cfg, data);
    cfg.trials = find((c==4) & prevc==5); % C-N
    dataout.data3R = ft_selectdata(cfg, data);
    
    cfg.trials = find((c==3) & (prevc==1 | prevc==4)); % IC-C
    dataout.data4L = ft_selectdata(cfg, data);
    cfg.trials = find((c==2) & (prevc==1 | prevc==4)); % IC-C
    dataout.data4R = ft_selectdata(cfg, data);
    cfg.trials = find((c==3) & (prevc==2 | prevc==3)); % IC-IC
    dataout.data5L = ft_selectdata(cfg, data);
    cfg.trials = find((c==2) & (prevc==2 | prevc==3)); % IC-IC
    dataout.data5R = ft_selectdata(cfg, data);
    cfg.trials = find((c==3) & prevc==5); % IC-N
    dataout.data6L = ft_selectdata(cfg, data);
    cfg.trials = find((c==2) & prevc==5); % IC-N
    dataout.data6R = ft_selectdata(cfg, data);
    
    cfg.trials = find(c==5 & (prevc==1 | prevc==4)); % N-C
    dataout.data7 = ft_selectdata(cfg, data);
    cfg.trials = find(c==5 & (prevc==2 | prevc==3)); % N-IC
    dataout.data8 = ft_selectdata(cfg, data);
    cfg.trials = find(c==5 & prevc==5); % N-N
    dataout.data9 = ft_selectdata(cfg, data);
    
%     fd = fieldnames(dataout);
%     for k=1:numel(fd)
%         dataout.(fd{k}).trialinfo(:, end+1) = k;
%     end
  otherwise
    error('unsupported reordering of data objects required');
end

