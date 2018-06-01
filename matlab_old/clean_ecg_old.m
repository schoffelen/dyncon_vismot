function [compecg] = clean_ecg(fname, trl, denoiseflag, list)

if nargin<3,
  denoiseflag = 0;
end

hdr   = read_header(fname);

%read in all data and do some filtering
cfg          = [];
cfg.trl      = trl;
cfg.datafile = fname;
if denoiseflag,
  %cfg.channel = channelselection({'MEG' '-A40' '-A107' '-A248' 'UACurrent'},hdr.label);
  cfg.channel = channelselection({'MEG' 'UACurrent'},hdr.label);
  cfg.denoise.channel    = channelselection('MEG', cfg.channel);
  cfg.denoise.refchannel = {'UACurrent'};
  cfg.denoise.hilbert = 'yes';
  data = preprocessing(cfg);
else
  cfg.channel  = {'MEG'};
end
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1;
if denoiseflag,
  if length(data.label)>size(data.trial{1},1),
    data.label = data.label(1:end-1); %remove UACurrent
  end
  cfg  = rmfield(cfg, 'datafile');
  cfg  = rmfield(cfg, 'trl');
  cfg  = rmfield(cfg, 'denoise');
  cfg.channel = {'MEG'};
  data = preprocessing(cfg, data);
else
  data = preprocessing(cfg);
end

megsel = match_str(data.label,channelselection('MEG', data.label));
if nargin<4,
  ecgsel = megsel;
else
  ecgsel = match_str(data.label,list);
end

for k = 1:length(data.trial)
  data.trial{k} = [data.trial{k}(megsel,:); mean(blc(data.trial{k}(ecgsel,:)))];
end
data.label = [data.label(megsel);'ECG'];

%create dummy ecg channel
ecg = zeros(1,max(trl(:)));
for k = 1:length(data.trial)
  ecg(trl(k,1):trl(k,2)) = data.trial{k}(end,:);
end
ecgorig  = ecg;
figure;plot(ecg(10000:50000));zoom
val      = input('peakindex?\n','s');
peakindx = str2num(val)+10000-1; %look it up yourself
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
