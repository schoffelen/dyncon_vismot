function [dataout,refdataout] = doPreprocessing(subject, condition, flag1, flag2, flag3)

if nargin==2,
  flag1 = 1; %collect meg-channels
  flag2 = 1; %collect ref-channels
  flag3 = 0; %do dftfiltering
elseif nargin==3,
  flag2 = 1;
  flag3 = 0;
elseif nargin==4,
  flag3 = 0; 
end

%switch denoising on if necessary
if isfield(subject, 'denoise') && subject.denoise==1, 
  denoiseflag = 1;
else
  denoiseflag = 0;
end

nrun = length(subject.runnames)
for rlop = 1:nrun
  %do semi-automatic artifact detection
  cd(subject.pathname);
  cd('artifact');
  load([subject.name,'eog-run',subject.runnames{rlop}(1)]);
  eog  = cfg.artfctdef.eog;
  load([subject.name,'jump-run',subject.runnames{rlop}(1)]);
  jump = cfg.artfctdef.jump;
  cd('../trl');
  load([subject.name,'trl-run',subject.runnames{rlop}(1),'cnd',num2str(condition)]);
  
  %add 4th column to trl matrix to code for condition
  cfg.trl(:, 4) = condition;

  %extract RT for the 6th column in the trl matrix
  startfix1 = cfg.trl(:,1);
  startfix2 = subject.startfix(:);
  [int,ia,ib] = intersect(startfix2, startfix1);
  cfg.trl(ib,6) = subject.rt(ia);

  %add 5th column for unique trial identifier
  cfg.trl(ib,5) = subject.trlid(ia);

  cfg.artfctdef.eog  = eog;
  cfg.artfctdef.jump = jump;
  
  %read header
  hdr = ft_read_header(cfg.datafile);

  %update trial specification and remove UACurrent and line noise
  cfg.channel         = {'MEG' '-A40' '-A107' '-A157' '-A248'};
  cfg.padding         = 5;
  cfg.dftfilter       = 'no';
  rayleigh            = hdr.Fs./round(5.*hdr.Fs);
  cfg.dftfreq         = round([49.8 50 50.2 99.8 100 100.2 150]./rayleigh).*rayleigh;
  cfg.detrend         = 'yes';
  cfg.artfctdef.reject = 'partial';
  cfg.artfctdef.minaccepttim = 0.25;
  cfg.artfctdef.type  = {'eog' 'jump'};
  cfg                 = ft_rejectartifact(cfg);
  if flag1,
    if flag3,
      cfg.dftfilter = 'yes';
    else
      cfg.dftfilter = 'no';
    end
    data                = ft_preprocessing(cfg);
    dataorig            = data;
  end
  if flag2,
    if flag3,
      cfg.dftfilter = 'yes';
    else
      cfg.dftfilter = 'no';
    end
    cfg.channel         = {'MEGREF'};
    refdata             = preprocessing(cfg);
  end
  cfg.dftfilter       = 'no';
  cfg.channel         = {'UACurrent'};
  UACdata             = preprocessing(cfg);
  if flag1, dataorig  = appenddata([], dataorig, UACdata); end
  if flag2, refdata   = appenddata([], refdata , UACdata); end
  if denoiseflag,
    cfg                 = [];
    cfg.denoise.refchannel = 'UACurrent';
    cfg.denoise.hilbert = 'yes';
    if flag1,
      cfg.denoise.channel = channelselection({'MEG'}, dataorig.label); 
      data                = preprocessing(cfg, dataorig);
    end
    if flag2,
      cfg.denoise.channel = channelselection({'MEGREF'}, refdata.label);
      refdata             = preprocessing(cfg, refdata); 
    end
  end

  %%do denoise_pca
  %dataorig       = data;
  %%cfg            = [];
  %%cfg.channel    = 'MEGREF';
  %%refdata        = preprocessing(cfg, dataorig);
  %%cfg.channel    = 'MEG';
  %%data           = preprocessing(cfg, dataorig);
  %%cfg.channel    = 'UACurrent';
  %%UACdata        = preprocessing(cfg, dataorig);
  %%clear dataorig; dataorig = data;
  %cfg            = [];
  %cfg.refchannel = {'MEGREFG'};
  %cfg.channel    = {'MEG'};
  %cfg.truncate   = 3;
  %data           = denoise_pca(cfg, dataorig, refdata);
  %clear dataorig; dataorig = data;
  %cfg            = [];
  %cfg.refchannel = {'MLzA' 'MLyA' 'MLxA'};
  %cfg.channel    = {'MEG'};
  %cfg.truncate   = 3;
  %data           = denoise_pca(cfg, dataorig, refdata);  
  %clear dataorig; dataorig = data;
  %cfg            = [];
  %cfg.refchannel = {'MLzaA' 'MLyaA' 'MLxaA'};
  %cfg.channel    = {'MEG'};
  %data           = denoise_pca(cfg, dataorig, refdata);
  
  %do rejectvisual
  %dataorig   = data;
  %cfg        = [];
  %cfg.method = 'summary';
  %cfg.channel = 'MEG';
  %cfg.keepchannel = 'yes';
  %data       = rejectvisual(cfg, dataorig);
  
  %if denoiseflag,
  %  trl1        = findcfg(data.cfg, 'trl');
  %  trl2        = findcfg(UACdata.cfg, 'trl');
  %  [int,i1,i2] = intersect(trl1, trl2, 'rows');
  %  UACdata.trial = UACdata.trial(i2);
  %  UACdata.time  = UACdata.time(i2);
  %  data          = appenddata([], data, UACdata);

  %  cfg                 = [];
  %  cfg.denoise.channel = channelselection({'MEG'}, dataorig.label); 
  %  cfg.denoise.refchannel = 'UACurrent';
  %  cfg.denoise.hilbert = 'yes';
  %  cfg.blc             = 'yes';
  %  data    = preprocessing(cfg, data);
  %
  %end
  if flag1, dataout{rlop}    = data;    end
  if flag2, refdataout{rlop} = refdata; end
end
if nrun==1,
  if flag1, dataout    = dataout{1};    end
  if flag2, refdataout = refdataout{1}; end
else
  if flag1, dataout    = appenddata([],dataout{:});    end
  if flag2, refdataout = appenddata([],refdataout{:}); end
end

if flag1 && ~flag2,
  refdataout = [];
end
if flag2 && ~flag1,
  dataout = [];
end


%cfg            = [];
%cfg.resamplefs = 500;
%cfg.detrend    = 'no';
%dataout        = resampledata(cfg, dataout);
