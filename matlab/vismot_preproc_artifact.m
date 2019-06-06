% This code was copied from doAnalysis and adapted to work with the current
% version of the pipeline.
subject = vismot_subjinfo(subjectname);

%% extract trial-definition for the respective conditions, and save
nconditions = 5;
for k = 1:length(subject.runnames)
 for m = 1:nconditions
   warning off;
   fname        = [subject.rawpath, '/', subject.name,'/', subject.scanname, subject.sessionname, subject.runnames{k},subject.datafile];
   cfg          = [];
   cfg.datafile = fname;
   cfg.trialfun = ['trialfun_condition',num2str(m)];
   cfg          = ft_definetrial(cfg);
   cd([subject.pathname,'trl/']);
   save([subject.name,'trl-run',subject.runnames{k}(1:end-1),'cnd',num2str(m)],'cfg');
 end
end

%% extract long trial-definition for detecting artifacts
for k = 1:length(subject.runnames)
 warning off;
 fname        = [subject.rawpath, '/', subject.name,'/', subject.scanname, subject.sessionname, subject.runnames{k},subject.datafile];
 cfg          = [];
 cfg.datafile = fname;
 cfg.trialfun = 'trialfun_longtrials';
 cfg          = ft_definetrial(cfg);
 cd([pathname,'trl/']);
 save([subjname,'trl-run',runnames{k}(1:end-1),'longtrials'],'cfg');
end


%% detect jumps
padding     = 5;
for k = 1:length(subject.runnames)
 warning off;
 load([subject.pathname,'/trl/',subject.name,'trl-run',subject.runnames{k}(1:end-1),'longtrials'],'cfg');
 trl   = cfg.trl;
 fname        = [subject.rawpath, '/', subject.name,'/', subject.scanname, subject.sessionname, subject.runnames{k},subject.datafile];
 cfg   = detect_jumps(fname, trl, padding, subject.denoise); 
 cd([subject.pathname,'artifact/']);
 save([subject.name,'jump-run',subject.runnames{k}(1:end-1)],'cfg');
end

%% compute topography of cardiogram IS THIS EVEN USED? HOW ARE ECG COMPONENTS REMOVED?

for k = 1:length(subject.runnames)
 %do semi-automatic artifact detection
 cd(subject.pathname);
 cd('artifact');
 load([subject.name,'jump-run',subject.runnames{k}(1)]);
 jump = cfg.artfctdef.jump;
 cd('../trl');
 load([subject.pathname,'trl/',subject.name,'trl-run',subject.runnames{k}(1:end-1),'longtrials'],'cfg');
 trl      = cfg.trl;
 dfile    = cfg.datafile; 
 cfg = [];
 cfg.artfctdef.type = 'jump';
 cfg.artfctdef.jump = jump;
 cfg.trl = trl;
 cfg.datafile = dfile;
 cfg  = ft_rejectartifact(cfg);

 trl   = cfg.trl;
  fname        = [subject.rawpath, '/', subject.name,'/', subject.scanname, subject.sessionname, subject.runnames{k},subject.datafile];
 load([pathname,'ecg/ecgtopolist']);
 %compecg = clean_ecg_old(fname, trl, SUBJ(j).denoise,list); 
 compecg = clean_ecg_old(fname, trl, subject.denoise); 
 cd([pathname,'ecg/']);
 save([subjname,'ecg-run',subject.runnames{k}(1:end-1)],'compecg');
end


%% detect eye-blinks
for k = 1:length(subject.runnames)
 warning off;
 load([subject.pathname,'/trl/',subject.name,'trl-run',subject.runnames{k}(1:end-1),'longtrials'],'cfg');
 try,
   load([subject.pathname,'/ecg/',subject.name,'ecg-run',subject.runnames{k}(1:end-1)],'compecg');
   ecg = compecg;
 catch
   ecg = [];
 end
 trl   = cfg.trl;
 fname = [subject.rawname,subject.runnames{k},subject.datafile];
 cfg   = detect_blinks(fname, trl, subject.denoise, ecg); 
 cd([pathname,'artifact/']);
 save([subjname,'eog-run',subject.runnames{k}(1:end-1)],'cfg');
end

