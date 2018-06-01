function [compecg] = clean_ecg2(data, list)

cfg         = [];
channel     = channelselection('MEG',data.label);
cfg.channel = channel;
cfg.numcomponent = length(cfg.channel);
cfg.method  = 'pca';
comp        = componentanalysis(cfg, data);

cfg        = [];
cfg.layout = '4D248.lay';
cfg.layout = prepare_layout(cfg);
componentbrowser(cfg, comp);
x = str2num(input('component number? ','s'));

cfg         = [];
cfg.channel = channel;
data        = preprocessing(cfg, data);
for k = 1:length(data.trial)
  data.trial{k} = [data.trial{k}; comp.trial{k}(x,:)];
end
data.label = [data.label;'ECG'];

%detect QRS-complexes
cfg                       = [];
cfg.artfctdef.ecg.inspect = {'A244';'A245';'A246'};
cfg                       = artifact_ecg(cfg, data);
trl                       = cfg.artfctdef.ecg.artifact;

if ~isempty(findcfg(data.cfg, 'resampletrl')),
  trlorig = findcfg(data.cfg, 'resampletrl');
else
  trlorig = findcfg(data.cfg, 'trl');
end

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
cfg.numcomponent = 100;
compecg        = componentanalysis(cfg, tlckecg);
