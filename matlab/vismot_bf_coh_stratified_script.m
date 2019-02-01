
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
subject = vismot_subjinfo(subjectname);
load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d10mm',subject.name)),'sourcemodel');

N = 75;
stratflag = 2;
[coh13,coh42,zx13,zx42,looptime] = vismot_bf_coh_stratified(subject,'sourcemodel',sourcemodel,'nrand',nrand,'stratflag', stratflag, 'N', N, 'lambda', lambda);

coh13 = rmfield(coh13, 'coh');
coh42 = rmfield(coh42, 'coh');

coh13.dcoh = coh13.coh_1-coh13.coh_2;
coh42.dcoh = coh42.coh_1-coh42.coh_2;

coh13 = rmfield(coh13, {'coh_1', 'coh_2'});
coh42 = rmfield(coh42, {'coh_1', 'coh_2'});

filename = fullfile(subject.pathname,'source', [subject.name, '_coh6d10mm_',sprintf('%03d', frequency)] );
save(filename, 'zx13', 'zx42', 'coh13', 'coh42');
