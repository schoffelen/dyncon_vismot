function [cfg] = detect_blinks(fname, trl, denoiseflag, ecg)

if nargin<3,
  denoiseflag = 0;
  ecg         = [];
elseif nargin<4,
  ecg         = [];
end

hdr   = ft_read_header(fname);
  
trllong      = trl;
trllong(:,1) = trllong(:,1)-round(0.5*hdr.Fs);
trllong(:,2) = trllong(:,2)+round(0.5*hdr.Fs);
trllong(:,3) = trl(:,3)-round(0.5*hdr.Fs);

%extract some channels which show blinks, and average
cfg          = [];
cfg.trl      = trllong;
cfg.datafile = fname;
cfg          = ft_checkconfig(cfg, 'dataset2files', 'yes');

if denoiseflag && isempty(ecg),
  fprintf('denoising the data\n');
  cfg.channel            = {'A174' 'A175' 'A176' 'A193' 'A194' 'A195' 'A227' 'A228' 'A247' 'A248' 'UACurrent'};
  cfg.denoise.channel    = {'A174' 'A175' 'A176' 'A193' 'A194' 'A195' 'A227' 'A228' 'A247' 'A248'};
  cfg.denoise.refchannel = {'UACurrent'};
  cfg.denoise.hilbert    = 'yes';
elseif denoiseflag && ~isempty(ecg),
  [a,b] = match_str(hdr.label, ecg.topolabel);
  cfg.channel            = [hdr.label(a);'UACurrent'];
  cfg.denoise.channel    = ecg.topolabel;
  cfg.denoise.refchannel = {'UACurrent'};
  cfg.denoise.hilbert    = 'yes';
  cfg.subspace           = eye(length(a)) - ecg.topo(b,1)*ecg.topo(b,1)';
  cfg.subspace(end+1,end+1) = 1;
elseif ~isempty(ecg),
  %do subspace projection
  [a,b] = match_str(hdr.label, ecg.topolabel);
  cfg.channel            = hdr.label(a);
  cfg.subspace           = eye(length(a)) - ecg.topo(b,1)*ecg.topo(b,1)';
else
  %no changes to cfg
end
data                   = ft_preprocessing(cfg);

cfg                          = [];
cfg.trl                      = trllong;
cfg.continuous               = 'yes';
cfg.artfctdef.eog.channel    = {'A174' 'A175' 'A176' 'A193' 'A194' 'A195' 'A227' 'A228' 'A247' 'A248'};
cfg.artfctdef.eog.boxcar     = 0.1;
cfg.artfctdef.eog.hilbert    = 'no';
cfg.artfctdef.eog.bpfilter   = 'yes';
cfg.artfctdef.eog.bpfreq     = [1 20];
cfg.artfctdef.eog.rectify    = 'yes';
cfg.artfctdef.eog.cutoff     = 4;
cfg.artfctdef.eog.fltpadding = 0;
cfg.artfctdef.eog.trlpadding = 0;
cfg.artfctdef.eog.artpadding = 0.1;
cfg.artfctdef.eog.feedback   = 'yes';
cfg.artfctdef.type           = 'eog';
cfg.artfctdef.reject         = 'partial';
cfg                          = ft_artifact_eog(cfg, data);
