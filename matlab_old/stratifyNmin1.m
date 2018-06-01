function [data] = stratifyNmin1(subject, flag)

%'randomly' select trials such that for each of the 4 response conditions
%the previous trials were matched for number (out of 5 total conditions) 
%of occurrences (and the previous had to be a correct trial)
%
%the output data contains n*4 trials for each of the response conditions.
%each of the response conditions consists of m*5 trials. each set designates
%the condition of the previous trial
%
%the order of the trials in the output array
%graphically: 1 2 3 4 5 6 ............          trial index
%             1 1 1 1 1 1 1 1 1 1 ....  4 4 4 4 response conditions
%             1 1 2 2 3 3 4 4 5 5 ....  4 4 5 5 condition of previous trial

if nargin==1,
  flag = 0;
end

cd(subject.pathname);
cd('data');
if flag==0,
  load([subject.name,'data']);
else
  load([subject.name,'data_aligned']);
end

cd(subject.pathname);
cd('stratifyRT');
load([subject.name,'stratifyRT']);

for k = 1:5
  warning off;
  if k==1,
    data = data1;
  elseif k==2,
    data = data2;
  elseif k==3,
    data = data3;
  elseif k==4,
    data = data4;
  elseif k==5,
    data = data5;
  end
  data = struct2double(data);
  warning on;
  trl  = findcfg(data.cfg, 'trl');
  if k<5,
    trl(:,4) = input{k}';
  else
    trl(:,4) = nan;
  end
  data.cfg.trl = trl;
  alldat{1,k}  = data;
end
data = ft_appenddata([], alldat{:});
clear alldat;

cfg        = [];
cfg.toilim = [-0.75 0-1./256];
cfg.minlength = 0.75-2/256;
data       = ft_redefinetrial(cfg, data);

cmat = subject.correct;
tmat = subject.trl;
tsmp = subject.trigsmp;
imat = zeros(size(tmat));
cond = zeros(size(tmat,1),3,5,5); 

tx = findcfg(data.cfg, 'trl');
tx = tx(:,1) - tx(:,3) + 1;
for k = 1:numel(tx)
  dd = abs(tx(k) - tsmp);
  imat(dd==min(dd(:))) = k;
end
%imat contains the index to the trial (in data) in the same
%format as how tmat is specified

%tmat contains for all trials in the experiment the condition
%row number corresponds to the 'true trial' (consisting of a
%sequence of 4 trials

for k = 1:5
  for m = 1:5
    cond(:,:,k,m) = tmat(:,2:4)==k & tmat(:,1:3)==m ...
                  & cmat(:,2:4)    & cmat(:,1:3)    ...
                  & imat(:,2:4)~=0;
  end %trial N-1
end %trial N
%cond contains for all 2nd to 4th trial in a row (210 times)
%a boolean 1 if
% -the present trial is of condition (3rd dim of cond) AND
% -the previous trial was correct AND
% -the previous trial was of condition (4th dim of cond)

%extract the minimum number of trials per cell
tmpimat = imat(:,2:end);
cond    = cond==1;
for k = 1:4
  for m = 1:5
    trlind{k,m} = tmpimat(cond(:,:,k,m));
  end
end

cellnum = cellfun(@numel, trlind);
minnum  = min(cellnum(:));

newtrlind = trlind;
for k = 1:numel(trlind)
  sel          = randperm(cellnum(k));
  newtrlind{k} = trlind{k}(sel(1:minnum));
end
selrpt = reshape(cell2mat(newtrlind)', [minnum*20 1]);

cfg         = [];
cfg.trials  = selrpt;
cfg.channel = 'MEG';
data        = ft_preprocessing(cfg, data);
