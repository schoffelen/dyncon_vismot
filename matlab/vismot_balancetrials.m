function dataout = vismot_balancetrials(datain)

% this function balances the data structures across rows and columns, based
% on the number of trials (and samples) of the second row of data
% structures, which are assumed to be the post stimulus data structures.
% the number of samples balancing is to account for RT differences, where
% the number of samples is taken as a proxy for the RT. This overall works,
% if not too many trials are cut short (or start late) due to partial
% artifacts. Also supports a single row of data structures.
%
% Note that for now the prestimulus data number of trials is only matched,
% not the number of samples (this is different than before).

if size(datain,1)>1
  indx = 2;
else
  indx = 1;
end

smp1    = cellfun('size',datain(indx,1).trial,2);
smp2    = cellfun('size',datain(indx,2).trial,2);
[indx1, indx2, nsmp1, nsmp2] = equalizeNsmp(smp1, smp2); % this function is in private

% adjust
tmpcfg = [];
tmpcfg.trials = indx1;
datain(indx,1)   = ft_selectdata(tmpcfg, datain(indx,1));

tmpcfg.trials = indx2;
datain(indx,2)   = ft_selectdata(tmpcfg, datain(indx,2));

tmpcfg = [];
tmpcfg.begsample = ones(numel(nsmp1),1);
tmpcfg.endsample = nsmp1(:);
dataout(indx,1)   = ft_redefinetrial(tmpcfg, datain(indx,1));
tmpcfg.begsample = ones(numel(nsmp2),1);
tmpcfg.endsample = nsmp2(:);
dataout(indx,2)   = ft_redefinetrial(tmpcfg, datain(indx,2));

% adjust the pre, if present
if size(datain,1)>1
  if numel(indx1)>=numel(datain(1,1).trial)
  else
    sel1 = sort(randperm(numel(indx1)));
    tmpcfg = [];
    tmpcfg.trials = sel1;
    datain(1,1) = ft_selectdata(tmpcfg, datain(1,1));
  end
  if numel(indx2)>=numel(datain(1,2).trial)
  else
    sel2 = sort(randperm(numel(indx2)));
    tmpcfg = [];
    tmpcfg.trials = sel2;
    datain(1,2) = ft_selectdata(tmpcfg, datain(1,2));
  end
  
  dataout(1,1) = datain(1,1);
  dataout(1,2) = datain(1,2);
end
