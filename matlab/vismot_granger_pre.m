function [granger, parcellation] = vismot_granger_pre(subject, varargin)

doprewhiten = istrue(ft_getopt(varargin, 'prewhiten',  false));
conditions  = ft_getopt(varargin, 'conditions', 1:5);
reverseflag = ft_getopt(varargin, 'reverseflag', 0);
split       = istrue(ft_getopt(varargin, 'split', false)); % compute individual conditions if true
label       = ft_getopt(varargin, 'label', 'all');

if ischar(subject)
	subject = vismot_subjinfo(subject);
end

[freqpre,tlckpre]      = vismot_spectral(subject,'output','fourier','conditions','previous','toi','pre');
if doprewhiten
  load(fullfile(subject.pathname,'data',[subject.name,'emptyroom']));
  
  cfgr        = [];
  cfgr.length = 0.5;
  noise       = ft_redefinetrial(cfgr, data);
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
  granger(k) = compute_granger(freq(k), parcellation, reverseflag);
end

function granger = compute_granger(freq, parcellation, reverseflag)

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
for k = 1:numel(begchunk)
	c1 = [1;1]*(begchunk(k):endchunk(k)).*2;
	c1(1,:) = c1(1,:)-1;
	n1 = size(c1,2);
	for m = k:numel(begchunk)
		c2 = [1;1]*(begchunk(m):endchunk(m)).*2;
	  c2(1,:) = c2(1,:)-1;
		n2 = size(c2,2);
	  
		% get the source level CSD matrix
		C = zeros(2*(n2+n1),2*(n2+n1),numel(freq.freq));
		f = F([c1(:);c2(:)],:);
		for p = 1:numel(freq.freq)
			C(:,:,p) = f*freq.crsspctrm(:,:,p)*f';
		end	
		
		for p = 1:n2
		  tic;
			fprintf('computing Granger for chunks %d/%d and %d/%d with reference parcel %d/%d\n',k,numel(begchunk)-1,m,numel(begchunk)-1,p,n2);
			cmbindx = [1:2:n1*2;2:2:n1*2]';
			cmbindx(:,3) = n1*2+p*2-1;
			cmbindx(:,4) = n1*2+p*2;
			tmpg = do_granger4x4(C, freq.freq, cmbindx);
			tmp  = cat(1,reshape(tmpg(1,2,:,:),[],numel(freq.freq)),reshape(tmpg(2,1,:,:),[],numel(freq.freq)));
			tmpcmb = parcellation.label(begchunk(k):endchunk(k));
			tmpcmb(:,2) = parcellation.label(begchunk(m)+p-1);
			tmpcmb = cat(1,tmpcmb,tmpcmb(:,[2 1]));
			
			for q = 1:size(tmpcmb,1)
				[~,indx] = match_str(tmpcmb(q,:), parcellation.label);
				G(indx(1),indx(2),:) = tmp(q,1:240);
			end
			toc;
		end		
	   
	end
end
warning on;

granger.grangerspctrm = G;
granger.freq          = freq.freq(1:240);
granger.label         = parcellation.label;
granger.dimord        = 'chan_chan_freq';
