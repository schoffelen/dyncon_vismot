function [event] = getEvent(subject, saveflag)

if nargin==1,
  saveflag = 1;
end

cd(subject.pathname);
cd('event');

if ~exist([subject.name,'event.mat'], 'file')
  for k = 1:length(subject.runnames)
    cd([subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname]);
    cd(subject.runnames{k});
    event{k} = read_event(subject.datafile);
  end
  
  if saveflag,
    cd(subject.pathname);
    cd('event');
    save([subject.name,'event'], 'event');
  end

else
  load([subject.name,'event']);
end
