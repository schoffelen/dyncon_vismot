function [oktrials,label,freqs] = stratifyFreqpre(freq)

label = freq{1}.grad.label;
cnt   = zeros(length(label),1);

cfg            = [];
cfg.channelcmb = {};
cfg.keeptrials = 'yes';
cfg.jackknife  = 'no';
for k = 1:length(freq)
  fd{k} = freqdescriptives(cfg, freq{k});  
  [i1,i2] = match_str(label, freq{k}.label);
  cnt(i1) = cnt(i1) + 1;
end
sel   = find(cnt==length(freq));
label = label(sel);

oktrials = cell(1,length(fd));
for k = 1:length(fd)
  siz         = size(fd{k}.powspctrm);
  oktrials{k} = logical(zeros([siz(1) length(label) siz(3:end)]));
end

for k = 1:length(label)
  fprintf('stratifying channel %s\n',label{k});
  for m = 1:length(fd{1}.freq)
    for kk = 1:length(fd)
      sel       = match_str(fd{kk}.label,label{k});
      input{kk} = log10(fd{kk}.powspctrm(:,sel,m))';
    end
    tmpcfg.equalbinavg = 'no';
    output = stratify(tmpcfg,input{:});
    for kk = 1:length(oktrials)
      oktrials{kk}(:,k,m) = isfinite(output{kk});
    end
  end
end
freqs = freq{1}.freq;
