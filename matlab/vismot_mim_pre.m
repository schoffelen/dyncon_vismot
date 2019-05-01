function [mim, parcellation] = vismot_mim_pre(subject, varargin)

doprewhiten   = istrue(ft_getopt(varargin, 'prewhiten',  false));
conditions    = ft_getopt(varargin, 'conditions', 1:5);
reverseflag   = ft_getopt(varargin, 'reverseflag', 0);
split         = istrue(ft_getopt(varargin, 'split', false)); % compute individual conditions if true
label         = ft_getopt(varargin, 'label', 'all');
dobaseline    = istrue(ft_getopt(varargin, 'dobaseline', false)); % only analyze trials that were not preceded by a previous trial, but by a baseline (only works in conditons previous)
doL1out       = istrue(ft_getopt(varargin, 'doL1out', false));
leaveouttrial = ft_getopt(varargin, 'leaveouttrial', false);

if doL1out && ~leaveouttrial
  error('please specify the trial number that should be left out.')
end

if ischar(subject)
	subject = vismot_subjinfo(subject);
end

[freqpre,tlckpre]      = vismot_spectral(subject,'output','fourier','conditions','previous','toi','pre', 'dobaseline', dobaseline, 'doL1out', doL1out, 'leaveouttrial', leaveouttrial);
if doprewhiten
  emptyroom = load(fullfile(subject.pathname,'data',[subject.name,'emptyroom']));
  
  cfgr        = [];
  cfgr.length = 0.5;
  noise       = ft_redefinetrial(cfgr, emptyroom.data);
  clear data;
  
  cfgd         = [];
  cfgd.detrend = 'yes';
  noise        = ft_preprocessing(cfgd, noise);
  
  cfg = [];
  cfg.covariance = 'yes';
  noise = ft_timelockanalysis(cfg, noise);

  for k = 1:numel(freqpre) 
    tmp1(k) = ft_denoise_prewhiten([], tlckpre(k), noise);
    tmp2(k) = ft_denoise_prewhiten([], freqpre(k), noise);
  end
  tlckpre = tmp1;
  freqpre = tmp2;
  clear tmp1 tmp2;
end

[source, parcellation] = vismot_bf_lcmv_pre(subject, tlckpre, 'truncate', 5);

for k = 1:numel(freqpre)
  tmp1(k) = ft_checkdata(freqpre(k), 'cmbrepresentation', 'fullfast');
end
freqpre = tmp1; clear tmp1;

if ~split
  % combine the crsspctrm
  freq           = freqpre(1);
  freq.crsspctrm = mean(cat(4,freqpre.crsspctrm),4);
else
  freq           = freqpre(conditions);
end

if ischar(label) && strcmp(label, 'all')
  % do nothing
else
  [a,b] = match_str(parcellation.label, label);
  parcellation.filter = parcellation.filter(a);
  parcellation.label = parcellation.label(a);
  parcellation.s = parcellation.s(a);
  parcellation.u = parcellation.u(a);
end


for k = 1:numel(freq)
  mim(k) = compute_mim(freq(k), parcellation, reverseflag);
end

function mim = compute_mim(freq, parcellation, reverseflag)

if reverseflag
	freq.crsspctrm = conj(freq.crsspctrm);
end

%parcellation.label = parcellation.label(1:10:end);

for k = 1:numel(parcellation.label)
	F((k-1)*2+(1:2),:) = parcellation.filter{k}(1:2,:);
end
 
n        = numel(parcellation.label);
chunk    = [(0:80:n) n];
begchunk = chunk(1:end-1)+1;
endchunk = chunk(2:end);

if 0
  % for testing
  begchunk = 1;
  endchunk = 50;
end

warning off;
G = zeros(n,n,240);

% compute the mim
M = zeros(size(F,1)/2,size(F,1)/2,240);
for p = 1:240
  tmp = F*freq.crsspctrm(:,:,p)*F';
  tmp = reshape(permute(reshape(tmp,[2 n 2 n]),[1 3 2 4]),[2 2 n^2]);
  tmp_r = real(tmp);
  tmp_i = imag(tmp);
  tmp_it = tmp_i;
  
  dum = tmp_i(1,2,:);
  tmp_it(2,1,:) = dum;
  dum = tmp_i(2,1,:);
  tmp_it(1,2,:) =dum;
 
  tmp_r = tmp_r(:,:,1:(n+1):end);
  tmp_rinv = inv2x2(tmp_r); % MvE: I have to cd to private in order for this to work.
  
  tmp_r1 = reshape(repmat(tmp_rinv,[1 1 1 n]),[2 2 n^2]);
  tmp_r2 = reshape(permute(repmat(tmp_rinv,[1 1 1 n]),[1 2 4 3]),[2 2 n^2]);
  
  cpath = pwd;
  cd('/project/3011085.03/scripts/fieldtrip/connectivity/private/')
  m = mtimes2x2(mtimes2x2(mtimes2x2(tmp_r1,tmp_i),tmp_r2),tmp_it);
  M(:,:,p) = reshape(m(1,1,:)+m(2,2,:),[n n]);
end
cd(cpath)

mim.mimspctrm     = M;
mim.freq          = freq.freq(1:240);
mim.label         = parcellation.label;
mim.dimord        = 'chan_chan_freq';
