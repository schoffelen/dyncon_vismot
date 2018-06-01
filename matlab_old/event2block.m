function [block] = event2block(subject)

%outputs the sample numbers (wrt raw data file)
%at which a new block commences in the second row
%first row corresponds with run-number

block = [];
for k = 1:length(subject.runnames)
  cd([subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname]);
  cd(subject.runnames{k});
  event = read_event(subject.datafile);
  type  = {event(:).type}';
  sel   = strmatch('TRIGGER', type);
  val   = [event(sel).value];
  smp   = [event(sel).sample];
  sel   = find(val==16); %begin of a trial
  smp   = smp(sel);
  dsmp  = diff([-smp(1) smp]);
  dsmp  = standardise(dsmp);
  sel   = find(dsmp>2);
  if sel(1) ~= 1,
    sel = [1 sel];
  end
  block = [block [k*ones(1,numel(sel)); smp(sel)]];
end
