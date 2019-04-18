addpath ~/matlab/toolboxes/export_fig/
load vismot_parcels
load groupresults_nmf; % in /project/3011085.03/analyse/granger

n = 94;


sel = 1;
nsel = size(c{sel},1);
indx = reshape(1:n^2,[],n);
indx(1:(n+1):end) = [];
S = zeros(nsel,n^2);
%S(:,indx) = c{sel};
S = c{sel};

A = a{sel};

g = [];
g.grangerspctrm = zeros(n);
g.grangerspctrm(mask(:)) = S(2,:);%reshape(S(1,:),[],n);
g.freq   = 0;
g.dimord = 'chan_chan';
g.label  = ulabel;
g.mask   = g.grangerspctrm;
srt = sort(g.grangerspctrm(:)); thr = srt(round(0.95.*(n^2)));
g.mask(g.mask<thr) = 0;
g.mask   = 2.*g.mask./max(g.mask(:));

figure('position',[100 100 1500 750],'color','w');
h1=subplot('position',[0.05 0.05 0.45 .9]);
h2=subplot('position',[0.55 0.1 0.4 .85]);
cfg.newfigure = 'no';

cmap = brewermap(64,'RdYlBu');
cmap = flipud(cmap);
cfg.colormap = cmap;

datadir = '/project/3011085.03/analysis/granger';
for sel = 1:numel(c)
  nsel = size(c{sel},1);
  indx = reshape(1:n^2,[],n);
  indx(1:(n+1):end) = []; 
  S = zeros(nsel,n^2);
  S(:,indx) = c{sel};
  %S = c{sel};
  A = a{sel};
  
  
  for k = 1:size(A,3)
    %figname = fullfile(datadir, sprintf('fig_nmf_ncomp%03d_%03d',sel+4, k));
  
    g.grangerspctrm = reshape(S(k,:),[],n);
    g.mask   = g.grangerspctrm;
    srt = sort(g.grangerspctrm(:)); thr = srt(round(0.8.*(n^2)));
    g.mask(g.mask<thr) = 0;

    g.mask   = 2.*g.mask./max(g.mask(:));
    
    axes(h1);hold on;ft_topoplotCC(cfg, g);
    axes(h1);ft_plot_layout(laynew, plotlayoptions{:});
    axis([-1.2 1.2 -1.2 1.2]);
    axes(h2);plot(0:0.5:119.5, squeeze(A(:,:,k))); abc = axis;
    hold on;plot(0:0.5:119.5, mean(mean(A(:,:,k),3),2), 'k', 'linewidth',2);
    axis([0 80 0 max(max(A(:,:,k)))]); xlabel('frequency (Hz)'); ylabel(sprintf('component %03d of %03d, with k = %d', k, size(A,3), sel+4));
    drawnow;
    %export_fig(figname, '-png');
    pause;
    cla(h1);cla(h2); clear abc
  end
end
