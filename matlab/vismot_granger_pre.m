function [granger, parcellation] = vismot_granger_pre(subject, reverseflag, split, conditions)

%explicitly add this path to circumvent some unresolved path issue when
%deploying through torque
%addpath('/opt/matlab/R2014b/toolbox/signal/signal');

if nargin < 2 || isempty(reverseflag),
	reverseflag = 0;
end

if nargin < 3 || isempty(split),
  split = 0;
end

if nargin < 4 || isempty(conditions),
  conditions = 1:5;
end

if ischar(subject),
	subject = vismot_subjinfo(subject);
end

[source, parcellation] = vismot_bf_lcmv_pre(subject);
if ~split
  [freqpre,tlckpre]      = vismot_spectral_pre(subject,'output','csd');
elseif split==1,
  % split data according to condition of previous trial, and compute
  % condition specific granger
  [freqpre,tlckpre]      =  vismot_spectral_pre(subject,'output','csd','conditions','previous');
elseif split>1,
  error('not supported');
end

if ~split
  % combine the crsspctrm
  freq           = freqpre(1);
  freq.crsspctrm = mean(cat(4,freqpre.crsspctrm),4);
else
  freq           = freqpre(conditions);
end

for k = 1:numel(freq)
  granger(k) = compute_granger(freq(k), parcellation, reverseflag);
end

function granger = compute_granger(freq, parcellation, reverseflag)

if reverseflag,
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

if 0,
  % for testing
  begchunk = 1;
  endchunk = 5;
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
			fprintf('computing Granger for chunks %d/%d and %d/%d with reference %d/%d\n',k,numel(begchunk),m,numel(begchunk),p,n2);
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
