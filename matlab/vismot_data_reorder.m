function [dataout] = vismot_data_reorder(datain, conditions)

switch conditions
  case 'current'
    % nothing needed
    dataout = datain;
  case 'previous'
    % order according to previous trial.
    
    % sanity check on input
    if numel(fieldnames(datain))~=5,
      error('unexpected number of data arguments in input');
    end
    
    data = ft_appenddata([],datain.data1,datain.data2,datain.data3,datain.data4,datain.data5);
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
    % add the condition of the previou trial to the trialinfo.
    data.trialinfo(:,end+1) = newc;
    
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
    error('not yet implemented');
  otherwise
    error('unsupported reordering of data objects required');
end

