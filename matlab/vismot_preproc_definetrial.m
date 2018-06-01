function [trl, trlclean] = vismot_preproc_definetrial(subject, condition)

nrun      = length(subject.runnames);
datafiles = vismot_subject2datafile(subject);
trl       = cell(nrun,1);
trlclean  = cell(nrun,1);
for rlop = 1:nrun
  %do semi-automatic artifact detection
  load(fullfile(subject.pathname,'artifact',[subject.name,'eog-run',subject.runnames{rlop}(1)]));
  eog  = cfg.artfctdef.eog;
  load(fullfile(subject.pathname,'artifact',[subject.name,'jump-run',subject.runnames{rlop}(1)]));
  jump = cfg.artfctdef.jump;
  load(fullfile(subject.pathname,'trl',[subject.name,'trl-run',subject.runnames{rlop}(1),'cnd',num2str(condition)]));
  trl{rlop} = cfg.trl;
  
  %add 4th column to trl matrix to code for condition
  trl{rlop}(:, 4) = condition;

  %extract RT for the 6th column in the trl matrix
  startfix1 = trl{rlop}(:,1);
  startfix2 = subject.startfix(:);
  [~,ia,ib] = intersect(startfix2, startfix1);
  trl{rlop}(ib,6) = subject.rt(ia);

  %add 5th column for unique trial identifier
  trl{rlop}(ib,5) = subject.trlid(ia);

  artfctdef.eog  = eog;
  artfctdef.jump = jump;
  
  cfg                  = [];
  cfg.headerfile       = datafiles{rlop};
  cfg.trl              = trl{rlop};
  cfg.artfctdef        = artfctdef;
  cfg.artfctdef.reject = 'partial';
  cfg.artfctdef.minaccepttim = 0.25;
  cfg.artfctdef.type  = {'eog' 'jump'};
  cfg                 = ft_rejectartifact(cfg);
  trlclean{rlop}      = cfg.trl;
end
