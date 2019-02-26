
if ~exist('frequency', 'var')
  error('frequency should be defined');
end

if ~exist('nrand', 'var')
  nrand = 100;
end

if ~exist('subjectname', 'var')
  error('subjectname needs to be defined');
end
% if ~exist('conditions', 'var')
%   error('conditions need to be specified (current, previous, or current_previous')
% end
if ~exist('smoothing', 'var')
  smoothing = 4;
end
if ~exist('lambda', 'var')
  lambda = '20%';
end
if ~exist('N', 'var')
  N = 75;
end
if ~exist('prewhiten', 'var')
  prewhiten = false;
end

subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d10mm',subject.name)),'sourcemodel');

if nrand>0
  % prune the inside, so that only the superficial dipoles are used
  dum = zeros(sourcemodel.dim);
  dum(sourcemodel.inside) = 1;
  dum = dum>0;
  se  = cat(3,[0 0 0;0 1 0;0 0 0],[0 1 0;1 1 1;0 1 0],[0 0 0;0 1 0;0 0 0])>0;
  dum2 = imerode(imerode(imerode(dum,se),se),se);
  sourcemodel.inside = (double(dum(:))-double(dum2(:)))>0;
end

% and then left/right symmetrize in anticipation of collapsing across
% left/right responses
dum = zeros(sourcemodel.dim);
dum(sourcemodel.inside) = 1;
dum2 = flip(dum,1);
dum = double(dum)+double(dum2);
sourcemodel.inside(:) = dum==2;


stratflag = 2;
[coh13,coh42,zx13,zx42,looptime] = vismot_bf_coh_stratified(subject,'sourcemodel',sourcemodel,'nrand',nrand,'stratflag', stratflag, 'N', N, 'lambda', lambda, 'prewhiten', prewhiten);

coh13 = rmfield(coh13, 'coh');
coh42 = rmfield(coh42, 'coh');

coh13.dcoh = coh13.coh_1-coh13.coh_2;
coh42.dcoh = coh42.coh_1-coh42.coh_2;

coh13 = rmfield(coh13, {'coh_1', 'coh_2'});
coh42 = rmfield(coh42, {'coh_1', 'coh_2'});

filename = fullfile(subject.pathname,'source', sprintf('%d/', N), [subject.name, '_coh6d10mm_',sprintf('%03d', frequency)] );
if istrue(prewhiten)
  filename = [filename '_prewhitened'];
end
save(filename, 'zx13', 'zx42', 'coh13', 'coh42', 'sourcemodel', 'N', 'nrand');
