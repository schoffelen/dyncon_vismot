function [tlck, covariance, dof] = tfr2tlck(freq, latency)

%convert tfr frequency structure with fourier into tlck structure with a covariance (real part of csd)
%for the time range latency, averaging over trials

if nargin==1,
  latency = 'all';
end

if isstr(latency) && strcmp(latency, 'all')
  latency = [freq.time(1) freq.time(end)];
end

ix                 = [nearest(freq.time, latency(1)):nearest(freq.time, latency(2))];
freq.fourierspctrm = permute(freq.fourierspctrm(:,:,:,ix), [1 4 2 3]); %make 'rpttap_time_chan_freq'
siz                = size(freq.fourierspctrm);
freq.fourierspctrm = reshape(freq.fourierspctrm, [siz(1)*siz(2) siz(3) siz(4)]);
siz                = size(freq.fourierspctrm);

covariance         = complex(zeros(siz(2),siz(2),siz(3)));
for k = 1:siz(3)
  tmp = transpose(freq.fourierspctrm(:,:,k));
  sel = find(isfinite(tmp(1,:)));
  tmp = tmp(:,sel);
  tmp = (tmp*tmp')./(numel(sel));
  covariance(:,:,k) = tmp;
  dof(k,1)          = numel(sel);
end

tlck       = [];
tlck.label = freq.label;
tlck.grad  = freq.grad;
tlck.cfg   = [];
tlck.time  = 0;
tlck.dimord = 'chan_time';
tlck.avg   = zeros(numel(tlck.label),1);
