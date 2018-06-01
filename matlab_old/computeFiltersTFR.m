function [source] = computeFiltersTFR(subject, frequency)

cd('/analyse/1/Project0002/tmpProject0030');
load([subject.name,'tfr2tlck']);
ix = nearest(freq, frequency);


allcov = covariance{1}(:,:,ix).*dof{1}(ix);
n      = dof{1}(ix);
for k = 2:4
  allcov = allcov + covariance{k}(:,:,ix).*dof{k}(ix);
  n      = n      + dof{k}(ix);
end
tlck.cov = real(allcov)./n;

%load vol
cd(subject.pathname);
cd('vol');
load([subject.name,'vol']);

%load grid and prune leadfields
[a,b] = match_str(tlck.label, tlck.grad.label);
load([subject.pathname,'grid/',subject.name,'grid6mm.mat']);
eval('newgrid = grid;');
for k = 1:length(newgrid.inside)
  indx  = newgrid.inside(k);
  tmplf = newgrid.leadfield{indx};
  newgrid.leadfield{indx} = tmplf(b,:);
end


cfg          = [];
cfg.method   = 'lcmv';
%cfg.subspace = 10;
cfg.lambda   = '10%';
cfg.vol      = vol;
cfg.grid     = newgrid;
cfg.fixedori = 'yes';
cfg.keepfilter = 'yes';
cfg.keepleadfield = 'yes';
cfg.feedback = 'textbar';
tlck.avg(:,2) = 0; %necessary to avoid error due to squeezing in sourceanalysis FIXME fix into reshape
tlck.time(2)  = 0.1; %irrelevant but keep consistent
source       = ft_sourceanalysis(cfg, tlck);
source       = sourcerestyle(source, 'new');
source       = estimate_fwhm_tetraNEW(source);

filt   = source.filter;
fwhm   = source.fwhm;
inside = source.inside;
