%% for the power results
alldir = '/project/3011085.03/';
if ~exist('stratifyflag', 'var'), stratifyflag = true; end
if ~exist('toi', 'var'), toi = 'post'; end
if ~exist('doroi', 'var'), doroi = false; end
if ~exist('dosave', 'var'), dosave = false; end

switch toi
  case 'post'
    if stratifyflag
      filename = fullfile([alldir, 'analysis/stat_bf_post_stratified.mat']);
    else
      filename = fullfile([alldir, 'analysis/stat_bf_post']);
    end
  case 'pre'
    if doroi
      filename = [alldir, 'analysis/stat_bf_pre.mat'];
    else
      filename = fullfile([alldir, 'analysis/stat_bf_pre_wholebrain']);
    end
end
load(filename, 'stat', 'c', 'ic')


if ~doroi
  % take the raw effect of the largest cluster
  for k=1:numel(stat)
    if isempty(stat{k}.posclusters), stat{k}.posclusters(1).prob = nan; end
    if isempty(stat{k}.negclusters), stat{k}.negclusters(1).prob = nan; end
    [~, ix] = nanmin([stat{k}.posclusters(1).prob stat{k}.negclusters(1).prob]);
    if ix==1
      tmpc(:,k) = mean(c(stat{k}.posclusterslabelmat==1,k,:),1);
      tmpic(:,k) = mean(ic(stat{k}.posclusterslabelmat==1,k,:),1);
    else
      tmpc(:,k) = mean(c(stat{k}.negclusterslabelmat==1,k,:),1);
      tmpic(:,k) = mean(ic(stat{k}.negclusterslabelmat==1,k,:),1);
    end
  end
  c = tmpc;
  ic = tmpic;
end

datC = [];
datC.time = 1:size(c,2);
datC.dimord = 'rpt_chan_time';
datC.label{1} = 'largest cluster';
datC.trial(:,1,:) = c;
datIC = removefields(datC, 'trial');
datIC.trial(:,1,:) = ic;


n=size(c,1);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_bayesfactor';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
bayesfactor = ft_timelockstatistics(cfg, datC, datIC);

if dosave
  save(filename, 'bayesfactor', '-append')
end

%% for the coherence results
alldir = '/project/3011085.03/';
filename = fullfile([alldir, 'analysis/stat_coh_pre.mat']);
load(filename, 'effect', 'idx_sign_uncorrected');

for k=1:numel(effect)
  c(:,k) = effect(k).rawcoh_C;
  ic(:,k) = effect(k).rawcoh_IC;
end

datC = [];
datC.time = 1:size(c,2);
datC.dimord = 'rpt_chan_time';
datC.label{1} = 'largest cluster';
datC.trial(:,1,:) = c;
datIC = removefields(datC, 'trial');
datIC.trial(:,1,:) = ic;


n=size(c,1);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_bayesfactor';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
bayesfactor = ft_timelockstatistics(cfg, datC, datIC);

for k=1:numel(effect)
  effect(k).bayesfactor = bayesfactor.bf10(k);
end
if dosave
  save(filename, 'effect', '-append');
end
effect(idx_sign_uncorrected).bayesfactor
