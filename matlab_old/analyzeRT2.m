function [rt,trl,correct,trigsmp,rs,runnr,trlid,startfix] = analyzeRT2(subject)

tmp   = getEvent(subject);
event = struct([]);
for k = 1:numel(tmp)
  if k==1,
    offset = 0;
  else
    offset = event(end).sample + 1000;
  end
  for m = 1:numel(tmp{k})
    tmp{k}(m).samplenew = tmp{k}(m).sample + offset;
    tmp{k}(m).run       = k;
  end
  if k==1,
    event = tmp{k};
  else
    event = cat(1, event, tmp{k});
  end
end


type   = {event(:).type}';
sel    = find(ismember(type, {'TRIGGER' 'RESPONSE'}));
event  = event(sel);
type   = type(sel);
value  = [event(:).value]';
value  = value - bitand(value,8192);
value  = value - bitand(value,4096)*(4095/4096);
sel    = find(value>0);
event  = event(sel);
value  = value(sel);
type   = type(sel);
sample = [event(:).sample]';
runix  = [event(:).run]';

%find trigger response sequences
sel      = [];
trlbegin = find(strcmp('TRIGGER', type(1:end-1)) & value(1:end-1)==8);
trlbegin = [trlbegin;length(type)];
for k = 1:length(trlbegin)-1
  tmpsel = trlbegin(k):trlbegin(k+1)-1;
  tmptr  = find(strcmp('TRIGGER', type(tmpsel)) & value(tmpsel)<=5);
  %some triggers are <x ms apart, which is probably an acquisition error
  tmpsel = tmptr-1+trlbegin(k);
  tmpsel(diff([0;sample(tmpsel)])<10) = [];
  if ~isempty(tmpsel),
    trl(k,:)  = value(tmpsel)';
    sel(k,:)  = tmpsel';
    trigsmp(k,:) = sample(tmpsel);
    startfix(k,:) = sample(tmpsel-1);
    runnr(k,:)   = runix(tmpsel);
    tmpresp   = find(strcmp('RESPONSE', type(tmpsel+1))); 
    resp(k,tmpresp) = value(tmpsel(tmpresp)+1)';
    rt(k,tmpresp)   = sample(tmpsel(tmpresp)+1)' - sample(tmpsel(tmpresp))';
  end
end
correct = (resp==32  & (trl==1 | trl==3)) | ...
          (resp==256 & (trl==2 | trl==4)) | ...
	  (resp==0   & trl==5);

ntrl = size(trl,1);
npt  = size(trl,2);

trlid = reshape((1:numel(trl)), [4 ntrl])';

%rt in trial N as a function of condition N
addpath /home/jan/matlab/robuststats
for k = 1:4
  %four conditions with responses
  sel        = intersect(find(trl==k),find(correct==1));
  rs.n(k)    = mean(rt(sel));
  rs.ntrim(k) = trimmean(rt(sel), 0.20);
  rs.ntrimsem(k)  = trimsem(rt(sel),  0.20);
end

%rt in trial N as a function of condition N-1
for k = 1:4
  for m = 1:5
    sel           = intersect(find(correct(:,2:end)==1),intersect(find(trl(:,2:end)==k),find(trl(:,1:end)==m)));
    tmprt         = rt(:,2:end);
    rs.nmin1(k,m) = mean(tmprt(sel));
    rs.nmin1trim(k,m) = trimmean(tmprt(sel),0.20);
    rs.nmin1trimsem(k,m)  = trimsem(tmprt(sel), 0.20);
    rs.sel{k,m}   = sel;
  end
end

%rt in the first one of a train of 'trials'
sel   = find(correct(:,1)==1 & trl(:,1)~=5);
tmprt = rt(sel, 1);
rs.nfirst        = mean(tmprt);
rs.nfirsttrim    = trimmean(tmprt, 0.2);
rs.nfirsttrimsem = trimsem(tmprt,  0.2); 





%delta curves
sel1   = intersect(find(trl==1),find(correct==1));
sel2   = intersect(find(trl==1+2), find(correct==1));
tmprt1 = rt(sel1); %congruent left
tmprt2 = rt(sel2); %incongruent left
srt1   = sort(tmprt1);
srt2   = sort(tmprt2);

%divide into quartiles
n1 = floor(linspace(0,numel(srt1),5));
n2 = floor(linspace(0,numel(srt2),5));
n1 = [n1(1:end-1)'+1 n1(2:end)'];
n2 = [n2(1:end-1)'+1 n2(2:end)'];
 
for k = 1:size(n1,1)
  ddat1(k,1) = mean([srt1(n1(k,1):n1(k,2));srt2(n2(k,1):n2(k,2))]);
  ddat1(k,2) = mean(srt1(n1(k,1):n1(k,2)))-mean(srt2(n2(k,1):n2(k,2)));
end

%delta curves
sel1   = intersect(find(trl==2+2),find(correct==1));
sel2   = intersect(find(trl==2), find(correct==1));
tmprt1 = rt(sel1); %congruent right
tmprt2 = rt(sel2); %incongruent right
srt1   = sort(tmprt1);
srt2   = sort(tmprt2);

%divide into quartiles
n1 = floor(linspace(0,numel(srt1),5));
n2 = floor(linspace(0,numel(srt2),5));
n1 = [n1(1:end-1)'+1 n1(2:end)'];
n2 = [n2(1:end-1)'+1 n2(2:end)'];
 
for k = 1:size(n1,1)
  ddat2(k,1) = mean([srt1(n1(k,1):n1(k,2));srt2(n2(k,1):n2(k,2))]);
  ddat2(k,2) = mean(srt1(n1(k,1):n1(k,2)))-mean(srt2(n2(k,1):n2(k,2)));
end

rs.ddatlft = ddat1;
rs.ddatrgt = ddat2;

%FIXME consider analysing results with respect to the fix-target interval
%according to Richard Ridderinkhof a strong pre cue with variable cue target
%intervals destroys the effect
