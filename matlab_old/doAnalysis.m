subjinfo;

%subjlist = 1;
subjlist = [1:3 5:6 8:16 18:20];
if 0,
  for j = subjlist
    subject = SUBJ(j);
    for condition = 1:4
      data = doPreprocessing(subject, condition, 1, 0 , 1);
      
      cfg  = [];
      cfg.toilim = [-0.4 -0.001];
      cfg.minlength = 0.38;
      data = ft_redefinetrial(cfg, data);
   
      %cfg  = [];
      %cfg.nneighbours = 50;
      %cfg.truncate    = 0.01;
      %cfg.channel     = {'MEG'};
      %data            = ft_denoise_sns(cfg, data);
   
      cfg  = [];
      cfg.blc = 'yes';
      cfg.detrend = 'no';
      cfg.resamplefs = 256;
      cfg.channel    = {'MEG'};
      data = ft_resampledata(cfg, data);
    
      warning off;
      if condition==1,
        data1 = struct2single(data);
      elseif condition==2,
        data2 = struct2single(data);
      elseif condition==3,
        data3 = struct2single(data);
      elseif condition==4,
        data4 = struct2single(data);
      end
      clear data;
    end
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/data/');
    save([SUBJ(j).name,'dataPre400'], 'data1', 'data2', 'data3', 'data4');
    clear data1 data2 data3 data4
  end
end

%subjlist = 1;
subjlist = [1:3 5:6 8:16 18:20];
if 0,
  for j = subjlist
    subject = SUBJ(j);
    for condition = 1:4
      data = doPreprocessing(subject, condition, 1, 0 , 1);
      
      cfg  = [];
      cfg.toilim = [0.25 0.75-1./data.fsample];
      cfg.minlength = 0.25;
      data = ft_redefinetrial(cfg, data);
   
      %cfg  = [];
      %cfg.nneighbours = 50;
      %cfg.truncate    = 0.01;
      %cfg.channel     = {'MEG'};
      %data            = ft_denoise_sns(cfg, data);
   
      cfg  = [];
      cfg.blc = 'yes';
      cfg.detrend = 'no';
      cfg.resamplefs = 256;
      cfg.channel    = {'MEG'};
      data = ft_resampledata(cfg, data);
    
      warning off;
      if condition==1,
        data1 = struct2single(data);
      elseif condition==2,
        data2 = struct2single(data);
      elseif condition==3,
        data3 = struct2single(data);
      elseif condition==4,
        data4 = struct2single(data);
      end
      clear data;
    end
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/data/');
    save([SUBJ(j).name,'dataPst500'], 'data1', 'data2', 'data3', 'data4');
    clear data1 data2 data3 data4
  end
end

if 1,
  for j = subjlist
    subject = SUBJ(j);
    freqpre = doFreqanalysisMtmfftPre400new(subject, 0, [2 20])  
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');
    save([subject.name,'freqPre400hanning'], 'freqpre');
    clear freqpre;
  end
end

if 0,
  for j = subjlist
    subject = SUBJ(j);
    freqpre = doFreqanalysisMtmfftPre400new(subject, 3.75, [8 40])  
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');
    save([subject.name,'freqPre400'], 'freqpre');
    clear freqpre;
  end
end

if 0,
  for j = subjlist
    subject = SUBJ(j);
    freqpre = doFreqanalysisMtmfftPre400new(subject, 7.5, [40 100])  
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');
    save([subject.name,'freqPre400high'], 'freqpre');
    clear freqpre;
  end
end

if 0,
  for j = subjlist
    subject = SUBJ(j);
    freqpst = doFreqanalysisMtmfftPst500new(subject, 4, [8 40])  
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');
    save([subject.name,'freqPst500low'], 'freqpst');
    clear freqpst;
  end
end

if 0,
  for j = subjlist
    subject = SUBJ(j);
    freqpst = doFreqanalysisMtmfftPst500new(subject, 8, [40 100])  
    cd('/data1/synchro1/Projects/JanMathijs/Project0030tmp/freq/');
    save([subject.name,'freqPst500'], 'freqpst');
    clear freqpst;
  end
end



%subjinfo;
%subjlist = [5:20];
%
%
%for j = subjlist
%
%subjname = SUBJ(j).name;
%rawname  = [SUBJ(j).rawpath,SUBJ(j).name,'/',SUBJ(j).scanname,SUBJ(j).sessionname];
%runnames = SUBJ(j).runnames;
%pathname = SUBJ(j).pathname;
%if ~isempty(SUBJ(j).datafile), 
%  datafile = SUBJ(j).datafile;
%else
%  datafile = 'c,rfDC';
%end
%
%if 0,
%%extract trial-definition for the respective conditions, and save
%nconditions = 5;
%for k = 1:length(runnames)
%  for m = 1:nconditions
%    warning off;
%    fname        = [rawname,runnames{k},datafile];
%    cfg          = [];
%    cfg.datafile = fname;
%    cfg.trialfun = ['trialfun_condition',num2str(m)];
%    cfg          = definetrial(cfg);
%    cd([pathname,'trl/']);
%    save([subjname,'trl-run',runnames{k}(1:end-1),'cnd',num2str(m)],'cfg');
%  end
%end
%end
%%------------------definetrial per condition
%
%if 0,
%%extract trial-definition for the long trials, and save
%%as of 07/04/2009 this has been re-run, because of an initial
%%bug in trialfun_longtrials. this means that the artifact 
%%detection has to be re-run as well
%for k = 1:length(runnames)
%  warning off;
%  fname        = [rawname,runnames{k},datafile];
%  cfg          = [];
%  cfg.datafile = fname;
%  cfg.trialfun = 'trialfun_longtrials';
%  cfg          = definetrial(cfg);
%  cd([pathname,'trl/']);
%  save([subjname,'trl-run',runnames{k}(1:end-1),'longtrials'],'cfg');
%end
%end
%%------------------definetrial long trials
%
%if 0,
%%detect jumps
%padding     = 5;
%for k = 1:length(runnames)
%  warning off;
%  load([pathname,'trl/',subjname,'trl-run',runnames{k}(1:end-1),'longtrials'],'cfg');
%  trl   = cfg.trl;
%  fname = [rawname,runnames{k},datafile];
%  cfg   = detect_jumps(fname, trl, padding, SUBJ(j).denoise); 
%  cd([pathname,'artifact/']);
%  save([subjname,'jump-run',runnames{k}(1:end-1)],'cfg');
%end
%end
%%---jumps
%
%if 0,
%%compute topography of cardiogram
%for k = 1:length(runnames)
%  %do semi-automatic artifact detection
%  cd(SUBJ(j).pathname);
%  cd('artifact');
%  load([SUBJ(j).name,'jump-run',SUBJ(j).runnames{k}(1)]);
%  jump = cfg.artfctdef.jump;
%  cd('../trl');
%  load([pathname,'trl/',subjname,'trl-run',runnames{k}(1:end-1),'longtrials'],'cfg');
%  trl      = cfg.trl;
%  dfile    = cfg.datafile; 
%  cfg = [];
%  cfg.artfctdef.type = 'jump';
%  cfg.artfctdef.jump = jump;
%  cfg.trl = trl;
%  cfg.datafile = dfile;
%  cfg  = rejectartifact(cfg);
%
%  trl   = cfg.trl;
%  fname = [rawname,runnames{k},datafile];
%  load([pathname,'ecg/ecgtopolist']);
%  %compecg = clean_ecg_old(fname, trl, SUBJ(j).denoise,list); 
%  compecg = clean_ecg_old(fname, trl, SUBJ(j).denoise); 
%  cd([pathname,'ecg/']);
%  save([subjname,'ecg-run',runnames{k}(1:end-1)],'compecg');
%end
%end
%%---ecg
%
%if 1,
%%detect eye-blinks
%for k = 1:length(runnames)
%  warning off;
%  load([pathname,'trl/',subjname,'trl-run',runnames{k}(1:end-1),'longtrials'],'cfg');
%  try,
%    load([pathname,'ecg/',subjname,'ecg-run',runnames{k}(1:end-1)],'compecg');
%    ecg = compecg;
%  catch
%    ecg = [];
%  end
%  trl   = cfg.trl;
%  fname = [rawname,runnames{k},datafile];
%  cfg   = detect_blinks(fname, trl, SUBJ(j).denoise, ecg); 
%  cd([pathname,'artifact/']);
%  save([subjname,'eog-run',runnames{k}(1:end-1)],'cfg');
%end
%end
%%---eyeblinks
%
%if 0,
%%preprocess data for timelocked analysis
%condition = 1;
%padding   = 5;
%for k = 1:length(runnames)
%  fname = [rawname,runnames{k},datafile];
%  hdr   = read_header(fname);
%  cd([pathname,'artifact/']);
%  load([subjname,'eog-run',runnames{k}(1:end-1)],'cfg');
%  eog  = cfg.artfctdef.eog;
%  %load([subjname,'jump-run',runnames{k}(1:end-1),'cnd',num2str(condition)],'cfg');
%  load([subjname,'jump-run',runnames{k}(1:end-1)],'cfg');
%  jump = cfg.artfctdef.jump;
%  
%  cd([pathname,'trl/']);
%  load([subjname,'trl-run',runnames{k}(1:end-1),'cnd',num2str(condition)],'cfg');
%  cfg.padding   = padding;
%  tmpdata       = preprocessing(cfg);
%  cfg.dftfilter = 'yes';
%  rayleigh      = hdr.Fs./round(padding*hdr.Fs);
%  dftbin        = round(50./rayleigh);
%  cfg.dftfreq   = [rayleigh*dftbin 2*rayleigh*dftbin];
%  cfg.dftinvert = 'yes';
%  cfg.channel   = 'MEG';
%  tmpdata2      = preprocessing(cfg);
%
%  cfg.artfctdef.jump = jump;
%  cfg.artfctdef.eog  = eog;
%  cfg.artfctdef.reject = 'partial';
%  cfg.artfctdef.minaccepttim = 0.25;
%  cfg.artfctdef.type = {'eog' 'jump'};
%  cfg.padding        = padding;
%  if k==1,
%    data  = rejectartifact(cfg, tmpdata);
%    data2 = rejectartifact(cfg, tmpdata2);
%  else
%    data  = appenddata([], data,  rejectartifact(cfg, tmpdata));
%    data2 = appenddata([], data2, rejectartifact(cfg, tmpdata2));
%  end
%  clear tmpdata tmpdata2;
%end
%cfg = [];
%cfg.truncate = 4;
%cfg.channel  = {'MEG' '-A40' '-A107' '-A157' '-A248'}; %known bad channels at the moment
%cfg.refchannel = {'MEG' '-A40' '-A107' '-A157' '-A248'}; %known bad channels at the moment
%data = denoise_pca(cfg, data, data2);
%clear data2;
%
%cfg        = [];
%cfg.method = 'summary';
%%cfg.metric = 'range';
%%data       = rejectvisual(cfg, data);
%cfg.metric = 'var';
%data       = rejectvisual(cfg, data);
%end
%
%end
