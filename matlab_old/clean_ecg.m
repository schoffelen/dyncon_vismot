function [compecg] = clean_ecg2(data, list)

if nargin==1,
  megsel = match_str(data.label,channelselection('MEG', data.label));
else
  megsel = match_str(data.label,list);
end


for k = 1:length(data.trial)
  data.trial{k} = [data.trial{k}(megsel,:); mean(blc(data.trial{k}(megsel,:)))];
end
data.label = [data.label(megsel);'ECG'];

%create dummy ecg channel
trl = findcfg(data.cfg, 'trl');
ecg = zeros(1,max(trl(:)));
for k = 1:length(data.trial)
  ecg(trl(k,1):trl(k,2)) = data.trial{k}(end,:);
end
ecgorig  = ecg;
figure;plot(ecg(20000:60000));zoom
val      = input('peakindex?\n','s');
peakindx = str2num(val)+20000-1; %look it up yourself
krn      = ecg(peakindx+[-500:500]);
ecg      = convn(ecg,fliplr(krn),'same');
for k = 1:length(data.trial)
  data.trial{k}(end,:) = ecg(trl(k,1):trl(k,2));
end

%detect QRS-complexes
cfg                       = [];
cfg.artfctdef.ecg.inspect = {'A244';'A245';'A246'};
cfg                       = artifact_ecg(cfg, data);
trl                       = cfg.artfctdef.ecg.artifact;

trlorig = findcfg(data.cfg, 'trl');
for k = 1:size(trl,1)
  tmp = find(trlorig(:,1)-trl(k,1) <=0 & trlorig(:,2)-trl(k,2)>=0)
  if ~isempty(tmp)
    ok(k) = 1;
  else
    ok(k) = 0;
  end
end

%redefinetrial
cfg          = [];
cfg.trl      = trl(find(ok),:);
cfg.trl(:,3) = 0;
dataecg      = redefinetrial(cfg, data);

%timelockanalysis
cfg          = [];
cfg.blc      = 'yes';
cfg.channel  = 'MEG';
tlckecg      = timelockanalysis(cfg, dataecg);

%componentanalysis
cfg            = [];
cfg.method     = 'pca';
cfg.numcomponent = 25;
compecg        = componentanalysis(cfg, tlckecg);
