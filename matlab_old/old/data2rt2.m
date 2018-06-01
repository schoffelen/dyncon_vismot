function [rt] = data2rt(data, event, flag)

%this function takes the trl matrix from a data structure
%and uses the corresponding event structure to return the
%reaction time per trial (so far implemented for freq data
%only. rt is expressed as number of samples

if nargin==2,
  flag = 1; %look in first column of trl for closest trigger
  %corresponding to cue onset
end

data = checkdata(data, 'datatype', 'freq');
trl  = findcfg(data.cfg, 'trl');
ntrl = size(trl,1);
if ntrl~=length(data.cumtapcnt), error('number of trl in trl ~= number of replicates in data'); end

if iscell(event)
  %cut trl into run-specific trl
  edges  = diff([1e10;trl(:,1)]);
  edges  = find(edges<0);
  edges(:,2) = [edges(2:end)-1;size(trl,1)];

  rt = [];
  for k = 1:length(event)
    tmp = trl2rt(trl(edges(k,1):edges(k,2),:),event{k}.event,flag);
    rt  = [rt;tmp];
  end
else
  rt = trl2rt(trl,event,flag);
end

function [rt] = trl2rt(trl,event,flag);

type = {event(:).type}';
trig = strmatch('TRIGGER',  type);
resp = strmatch('RESPONSE', type);
sel  = sort([trig;resp]);
type = type(sel);
trig = strmatch('TRIGGER',  type);
resp = strmatch('RESPONSE', type);
val  = [event(sel).value]';
smp  = [event(sel).sample]';

trigval = val(trig);
trigsmp = smp(trig);
respval = val(resp);
respsmp = smp(resp);

trigval(find(bitand(trigval,4096))) = trigval(find(bitand(trigval,4096)))-4095;
trigval(find(trigval<=4))           = 1; %convert all cue-onsets to 1

cuesmp  = trigsmp(trigval==1);

ntrl = size(trl,1);
for k = 1:ntrl
  if flag==1,
    %cue = cuesmp(nearest(cuesmp,trl(k,flag)));
    cue = cuesmp(find(cuesmp-trl(k,flag)<0,1,'last')+1);
  elseif flag==2,
    cue = cuesmp(find(cuesmp-trl(k,flag)>0,1,'first'));
  end
  press   = respsmp(find(respsmp>cue,1));
  rt(k,1) = press - cue;
end
