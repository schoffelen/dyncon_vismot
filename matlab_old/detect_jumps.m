function [cfg] = detect_jumps(fname, trl, padding, denoiseflag)

if nargin<4, 
  denoiseflag = 0;
elseif nargin<3,
  denoiseflag = 0;
  padding     = 5;
end

hdr   = ft_read_header(fname);
  
trllong      = trl;
trllong(:,1) = trllong(:,1)-round(0.5*padding*hdr.Fs);
trllong(:,2) = trllong(:,2)+round(0.5*padding*hdr.Fs);
trllong(:,3) = trl(:,3)-round(0.5*padding*hdr.Fs);

cfg          = [];
cfg.datafile = fname;
cfg          = dataset2files(cfg);
cfg.trl      = trllong;
%cfg.trl      = trl;
cfg.blc      = 'yes';
if denoiseflag,
  cfg.channel            = {'MEG' 'MEGREF' 'UACurrent'};
  cfg.denoise.channel    = ft_channelselection({'MEG' 'MEGREF'},hdr.label);
  cfg.denoise.refchannel = {'UACurrent'};
  cfg.denoise.hilbert    = 'yes';
else
  cfg.channel            = {'MEG' 'MEGREF'};
end
%data         = preprocessing(cfg);
%data         = struct2single(data);
%if denoiseflag,
%  data.label = data.label(1:end-1); 
%end
%cfg                   = [];
cfg.artfctdef         = [];
cfg.artfctdef.type    = 'jump';
cfg.artfctdef.reject  = 'partial';
cfg.artfctdef.jump.artpadding = 2.5;
cfg.artfctdef.jump.cutoff     = 25;
cfg.artfctdef.jump.feedback   = 'yes';
if denoiseflag,
  cfg.artfctdef.jump.denoise = cfg.denoise;
  cfg.artfctdef.jump.channel = {'MEG' 'MEGREF' 'UACurrent'};
else
  cfg.artfctdef.jump.channel = {'MEG' 'MEGREF'};
end
cfg.trl               = trl;
cfg.continuous        = 'yes';
%cfg                   = artifact_jump(cfg, data);
cfg                   = ft_artifact_jump(cfg);

%keyboard
%cfg          = [];
%cfg.method   = 'summary';
%cfg.metric   = 'range';
%cfg.channel  = 'MEG';
%cfg.medianfilter  = 'yes';
%cfg.medianfiltord = 9;
%cfg.absdiff  = 'yes';
%data         = rejectvisual(cfg, data);
%
%cfg          = struct2double(data.cfg);
%cfg.trl(:,1) = cfg.trl(:,1)+round(0.5*padding*hdr.Fs);
%cfg.trl(:,2) = cfg.trl(:,2)-round(0.5*padding*hdr.Fs);
%cfg.trl(:,3) = cfg.trl(:,3)+round(0.5*padding*hdr.Fs);
