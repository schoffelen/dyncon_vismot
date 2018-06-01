function [something] = processEmptyroom(subject)

cd(subject.rawpath);
cd(subject.name);
cd(subject.scanname);
cd(subject.sessionname);
cd(subject.emptyrun{1});

hdr = read_header(subject.datafile);

cfg                      = [];
cfg.trialfun             = 'trialfun_general';
cfg.trialdef.triallength = 0.5;
cfg.datafile             = subject.datafile;
cfg.dataset              = subject.datafile;
cfg                      = definetrial(cfg);

cfg.channel = {'MEG' 'MEGREF'};
data        = preprocessing(cfg);

cd(subject.pathname);
cd('data');
load([subject.name,'data'],'data1');
weights = findcfg(data1.cfg, 'pca');

cfg     = [];
cfg.pca = weights;
data2   = denoise_pca(cfg, data);

%cfg         = [];
%cfg.channel = 'MEG';
%data        = preprocessing(cfg, data);

cfg            = [];
cfg.resamplefs = 256;
cfg.blc        = 'yes';
cfg.detrend    = 'no';
data2          = resampledata(cfg, data2);
for k = 1:length(data2.trial)
  %make 128 samples long
  data2.trial{k} = data2.trial{k}(:,1:end-1);
  data2.time{k}  = data2.time{k}(1:end-1);
end

cfg        = [];
cfg.method = 'summary';
data2      = rejectvisual(cfg, data2);

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [6 36];
cfg.pad    = 128./data2.fsample;
freq       = freqanalysis(cfg, data2);

datcov = zeros(248,248,length(freq.freq));
for k = 1:length(freq.freq)
  freq.freq(k)
  tmp = selectdata(freq,'foilim',freq.freq(k));
  tmp = checkdata(tmp, 'cmbrepresentation', 'full');
  datcov(:,:,k) = squeeze(mean(tmp.crsspctrm));
end
