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
  addpath ~/matlab/toolboxes/export_fig/
  
  load gaTFR_high
  
  cfg = [];
  cfg.frequency = [38 100];
  cfg.latency   = [-0.4 0.5];
  Favg_high = ft_selectdata(cfg, Favg);

  load gaTFR_low
  Favg_low = Favg;
  
  cfg = [];
  cfg.colormap = '*RdBu'; % this requires the currently open PR of JM's colormap branch
  cfg.zlim = 'maxabs';
  cfg.xlim = [0.2 inf];
  cfg.layout = '4D248_helmet.mat';
  cfg.gridscale = 150;
  figure;ft_topoplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_topo_high_allfreq');
  cfg.ylim = [34 50];
  cfg.highlightchannel = {'A140' 'A167' 'A168' 'A189'};
  cfg.highlightcolor = [0 0 0];
  cfg.highlightsize = 24;
  cfg.highlightsymbol = '.';
  cfg.highlight = 'on';
  figure;ft_topoplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_topo_high_lowgamma');
  cfg.ylim = [54 66];
  cfg.highlightchannel = {'A134' 'A161' 'A162' 'A183'};
  cfg.highlightcolor = [0 0 0];
  figure;ft_topoplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_topo_high_midgamma');
  cfg.ylim = [74 86];
  cfg.highlightchannel = {'A8' 'A9' 'A23' 'A24'};
  cfg.highlightcolor = [1 1 1];
  figure;ft_topoplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_topo_high_highgamma');
  cfg.ylim = [5.5 6.5];
  cfg.xlim = [0.1 inf];
  cfg.highlightchannel = {'A24' 'A25' 'A43' 'A44'};
  cfg.highlightcolor = [0 0 0];
  figure;ft_topoplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_topo_low_theta');
  cfg.ylim = [9.5 10.5];
  cfg.highlightchannel = {'A188' 'A205' 'A206' 'A222'};
  cfg.highlightcolor = [1 1 1];
  figure;ft_topoplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_topo_low_alpha');
  cfg.ylim = [21.5 22.5];
  cfg.highlightchannel = {'A25' 'A26' 'A44' 'A45'};
  cfg.highlightcolor = [1 1 1];
  figure;ft_topoplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_topo_low_beta');
  
  cfg = rmfield(cfg, 'xlim');
  cfg = rmfield(cfg, 'ylim');
  
  cfg.channel = {'A140' 'A167' 'A168' 'A189'};
  cfg.title   = 'right occipito-parietal';
  figure;ft_singleplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_lowgamma');
  
  cfg.channel = {'A134' 'A161' 'A162' 'A183'};
  cfg.title   = 'left occipito-parietal';
  figure;ft_singleplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_midgamma');
  
  cfg.channel = {'A8' 'A9' 'A23' 'A24'};
  cfg.title   = 'central';
  figure;ft_singleplotTFR(cfg, Favg_high);
  export_fig('-png','-eps','m2', 'tfr_highgamma');
  
  cfg.channel = {'A24' 'A25' 'A43' 'A44'};
  cfg.title   = 'central';
  cfg.zlim = [-0.4 0.4];
  figure;ft_singleplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_theta');
  
  cfg.channel = {'A188' 'A205' 'A206' 'A222'};
  cfg.title   = 'occipito-parietal';
  figure;ft_singleplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_alpha');
  
  cfg.channel = {'A25' 'A26' 'A44' 'A45'};
  cfg.title   = 'sensorimotor';
  figure;ft_singleplotTFR(cfg, Favg_low);
  export_fig('-png','-eps','m2', 'tfr_beta');
  
  
end


