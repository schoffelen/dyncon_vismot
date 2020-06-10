datadir = '/project/3011085.03/analysis/freq';
cd(datadir);

docompute = false;
if docompute
  d = dir('*high.mat');
  for k = 1:numel(d)
    load(d(k).name);
    F{k} = freq(1);
    F{k}.powspctrm = squeeze(nanmean(cat(1,freq(:).powspctrm),1));
    F{k}.dimord = 'chan_freq_time';
    
    f13{k} = freq(1);
    f13{k}.powspctrm = squeeze(nanmedian(freq(1).powspctrm)./nanmedian(freq(3).powspctrm))-1;
    f13{k}.dimord = 'chan_freq_time';
    
    f42{k} = freq(4);
    f42{k}.powspctrm = squeeze(nanmedian(freq(4).powspctrm)./nanmedian(freq(2).powspctrm))-1;
    f42{k}.dimord = 'chan_freq_time';
    
  end
  
  cfgb = [];
  cfgb.baselinetype = 'relchange';
  cfgb.baseline = [-0.5 -0.2];
  for k = 1:numel(F)
    Fb{k} = ft_freqbaseline(cfgb, F{k});
  end
  
  cfg=[];
  cfg.appenddim = 'rpt';
  cfg.parameter = 'powspctrm';
  Favg = ft_appendfreq(cfg, Fb{:});
  f13avg = ft_appendfreq(cfg, f13{:});
  f42avg = ft_appendfreq(cfg, f42{:});
  
  
  cfg = [];
  cfg.avgoverrpt = true;
  Favg = ft_selectdata(cfg, Favg);
  f13avg = ft_selectdata(cfg, f13avg);
  f42avg = ft_selectdata(cfg, f42avg);
  
  cfgp = [];
  cfgp.zlim = 'maxabs';
  cfgp.layout = '4D248_helmet.mat';
  %figure;ft_topoplotTFR(cfgp, Favg);
  
  
  load labels_mirrored
  
  sel = intersect(match_str(label(:,1),f13avg.label),match_str(label(:,2),f13avg.label));
  label = label(sel,:);
  [~,reordered] = match_str(label(:,1),label(:,2));
  
  sel = match_str(f13avg.label, unique(label(:)));
  
  
  favg = f13avg;
  favg.powspctrm = (f13avg.powspctrm(sel,:,:) + f42avg.powspctrm(reordered,:,:))./2;
  favg.label = f13avg.label(sel);
  save('gaTFR_high', 'Favg', 'f13avg', 'f42avg', 'favg');
end


doplot = true;
if doplot
  load gaTFR_high
  
  cfg = [];
  cfg.frequency = [34 100];
  cfg.latency   = [-0.4 0.5];
  Favg_high = ft_selectdata(cfg, Favg);

  load gaTFR_low
  Favg_low = Favg;
  
  cfg = [];
  cfg.colormap = '*RdBu'; % this requires the currently open PR of JM's colormap branch
  cfg.zlim = 'maxabs';
  cfg.layout = '4D248_helmet.mat';
  cfg.gridscale = 150;
  ft_topoplotTFR(cfg, Favg_high);
  
  
end


