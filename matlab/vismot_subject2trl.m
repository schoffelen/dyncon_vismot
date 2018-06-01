function trl = vismot_subject2trl(subject,type)

if nargin<2
  type = 'cnd';
end

if ischar(subject)
  subject = vismot_subjinfo(subject);
end

trl = cell(1, numel(subject.runnames));
for k = 1:numel(trl)
  trl{k} = zeros(0,5);
end
runnumbers = zeros(1, numel(trl));
for k = 1:numel(trl)
  runnumbers(1,k) = str2double(subject.runnames{k}(1));
end

switch type
  case {'cnd' 'longtrials'}
    trldir = fullfile(subject.pathname, 'trl');
    d      = dir([trldir,filesep,subject.name,'*',type,'*']);
    for k = 1:numel(d)
      if strcmp(type, 'cnd')
        load(fullfile(trldir,d(k).name));
        
        runnr = str2double(d(k).name(strfind(d(k).name,'run')+3));
        cndnr = str2double(d(k).name(strfind(d(k).name,'cnd')+3));
        
        tmp = cfg.trl;
        tmp(:,end+1) = cndnr; % append this info to the trl matrix
        tmp(:,end+1) = runnr;
        
        indx = find(runnumbers==runnr);
        trl{indx} = cat(1,trl{indx}, tmp);   
      else
        error('not yet implemented');
      end
    end
  otherwise
    error('unknown trial type requested');
end
