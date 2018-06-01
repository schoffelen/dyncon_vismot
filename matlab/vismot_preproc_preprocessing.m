function [dataout] = vismot_preproc_preprocessing(subject, trl)

% preprocess the data and do a proper dftfilter, use trl specification in trl

if ischar(subject)
  subject = vismot_subjinfo(subject);
end

% switch denoising on if necessary
if isfield(subject, 'denoise') && subject.denoise==1 
  denoiseflag = 1;
else
  denoiseflag = 0;
end

nrun       = numel(subject.runnames);
dataout    = cell(nrun,1);
refdataout = cell(nrun,1);

%convert trl to cell array of matrices, consistent with the runs
if nrun>1
  run      = [0;find(diff(trl(:,1))<0)]+1;
  run(:,2) = [run(2:end,1)-1;size(trl,1)];
  trlold   = trl;
  clear trl;
  for m = 1:nrun
    trl{m,1} = trlold(run(m,1):run(m,2),:);
  end
else
  if ~iscell(trl), trl = {trl}; end;
end

datafiles = vismot_subject2datafile(subject);

for rlop = 1:nrun
  % read header
  cfg          = [];
  cfg.datafile = datafiles{rlop};
  hdr          = ft_read_header(cfg.datafile);

  % update trial specification and remove UACurrent and line noise
  cfg.trl       = trl{rlop};
  cfg.channel   = {'MEG' '-A40' '-A107' '-A157' '-A248'};
  cfg.padding   = 5;
  cfg.dftfilter = 'yes';
  cfg.demean    = 'yes';
  rayleigh      = hdr.Fs./round(5.*hdr.Fs);
  cfg.dftfreq   = round([49.8 50 50.2 99.8 100 100.2 148.8 150 150.2]./rayleigh).*rayleigh;
  data          = ft_preprocessing(cfg);
  dataorig      = data;
  cfg.channel   = {'MEGREF'};
  refdata       = ft_preprocessing(cfg);
  cfg.dftfilter = 'no';
  cfg.channel   = {'UACurrent'};
  UACdata       = ft_preprocessing(cfg);
  dataorig      = ft_appenddata([], dataorig, UACdata);
  refdata       = ft_appenddata([], refdata , UACdata);
  
  % this is for the removal of the coil-signal
  if denoiseflag
    cfg                 = [];
    cfg.denoise.refchannel = 'UACurrent';
    cfg.denoise.hilbert = 'yes';
    cfg.denoise.channel = ft_channelselection({'MEG'}, dataorig.label); 
    data                = ft_preprocessing(cfg, dataorig);
    cfg.denoise.channel = ft_channelselection({'MEGREF'}, refdata.label);
    refdata             = ft_preprocessing(cfg, refdata); 
  end

  dataout{rlop}    = data;
  refdataout{rlop} = refdata; 
end
if nrun==1
  dataout    = dataout{1};
  refdataout = refdataout{1};
else
  % average the gradiometer definition
  for k = 1:nrun
    if isfield(dataout{k}, 'hdr')
      w(k,1)    = dataout{k}.hdr.nSamples;
    else
      w(k,1)    = 1;
    end
    sens(k,1) = dataout{k}.grad;
  end
  if any(w==1)
    w(:) = 1;
  end
  grad = ft_average_sens(sens, 'weights', w);
  
  dataout    = ft_appenddata([],dataout{:});   
  refdataout = ft_appenddata([],refdataout{:});

  dataout.grad    = grad;
  refdataout.grad = grad;
end

cfg          = [];
cfg.truncate = 4;
dataout      = ft_denoise_pca(cfg, dataout, refdataout);

cfg          = [];
cfg.resamplefs = 300;
dataout      = ft_resampledata(cfg, dataout);
