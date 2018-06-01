function [dataout,refdataout] = doPreprocessing(subject, trl)

%re-preprocess the data and do a proper dftfilter, use trl specification in trl

%switch denoising on if necessary
if isfield(subject, 'denoise') && subject.denoise==1, 
  denoiseflag = 1;
else
  denoiseflag = 0;
end

nrun = length(subject.runnames);
%FIXME do something here with the trl matrix
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

for rlop = 1:nrun
  %if nrun>1, error('nrun>1 not yet supported'); end
  %read header
  cfg = [];
  cfg.datafile = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{rlop},subject.datafile];
  hdr = read_header(cfg.datafile);

  %update trial specification and remove UACurrent and line noise
  cfg.trl             = trl{rlop};
  cfg.channel         = {'MEG'};
  cfg.padding         = 5;
  cfg.dftfilter       = 'yes';
  rayleigh            = hdr.Fs./round(5.*hdr.Fs);
  cfg.dftfreq         = round([49.8 50 50.2 99.8 100 100.2 150]./rayleigh).*rayleigh;
  cfg.detrend         = 'yes';
  data                = preprocessing(cfg);
  dataorig            = data;
  cfg.channel         = {'MEGREF'};
  refdata             = preprocessing(cfg);
  cfg.dftfilter       = 'no';
  cfg.channel         = {'UACurrent'};
  UACdata             = preprocessing(cfg);
  dataorig  = appenddata([], dataorig, UACdata);
  refdata   = appenddata([], refdata , UACdata);
  if denoiseflag,
    cfg                 = [];
    cfg.denoise.refchannel = 'UACurrent';
    cfg.denoise.hilbert = 'yes';
    cfg.denoise.channel = channelselection({'MEG'}, dataorig.label); 
    data                = preprocessing(cfg, dataorig);
    cfg.denoise.channel = channelselection({'MEGREF'}, refdata.label);
    refdata             = preprocessing(cfg, refdata); 
  end

  dataout{rlop}    = data;
  refdataout{rlop} = refdata; 
end
if nrun==1,
  dataout    = dataout{1};
  refdataout = refdataout{1};
else
  dataout    = appenddata([],dataout{:});   
  refdataout = appenddata([],refdataout{:});
end
