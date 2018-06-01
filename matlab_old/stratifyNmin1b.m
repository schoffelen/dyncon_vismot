function [output] = selectTrials(subject, type, fname)

%'randomly' select trials such that for each of the 4 response conditions
%the previous trials were matched for number (out of 5 total conditions) 
%of occurrences (and the previous had to be a correct trial)


%%
%%the output data contains n*4 trials for each of the response conditions.
%%each of the response conditions consists of m*5 trials. each set designates
%the condition of the previous trial
%%
%%the order of the trials in the output array
%%graphically: 1 2 3 4 5 6 ............          trial index
%%             1 1 1 1 1 1 1 1 1 1 ....  4 4 4 4 response conditions
%%             1 1 2 2 3 3 4 4 5 5 ....  4 4 5 5 condition of previous trial

if nargin<3,
  fname = [subject.pathname,'freq/',subject.name,'mtmfftPre000'];
end

if nargin<2,
  error('type of stratification is required in input');
end

load(fname);
for k = 1:4
  trl{1,k} = findcfg(freqpre{k}.cfg, 'trl');
end
newtrl = fixTrl(subject, trl); clear trl;
trl    = newtrl; clear newtrl;

cmat = subject.correct;
tmat = subject.trl;
tsmp = subject.trigsmp;
imat = zeros(size(tmat));
cond = zeros(size(tmat,1),3,5,5); 

tx = trl;
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

%cellnum = cellfun(@numel, trlind);
%minnum  = min(cellnum(:));
%
%newtrlind = trlind;
%for k = 1:numel(trlind)
%  sel          = randperm(cellnum(k));
%  newtrlind{k} = trlind{k}(sel(1:minnum));
%end
%selrpt = reshape(cell2mat(newtrlind)', [minnum*20 1]);

if iscell(type)
  output = seltrl(cmat, imat, tmat, type{1});
  for k = 2:numel(type)
    output = intersect(output, seltrl(cmat, imat, tmat, type{k}));
  end
else
  output = seltrl(cmat, imat, tmat, type);
end

function [output] = seltrl(cmat, imat, tmat, type)

switch type
  case 'previouscorrect'
    tmpcmat = [ones(size(cmat,1),1) cmat(:,1:3)];
    output  = sort(imat(logical(tmpcmat))); 
    output  = output(output>0);
  case 'previouscongruent'
    tmptmat(:,2:4) = tmat(:,1:3);
    tmptmat(:,1)   = 0;
    %[1 4] are congruent
    tmptmat = tmptmat==1 | tmptmat==4;
    output  = sort(imat(tmptmat));
    output  = output(output>0);
  case 'previousincongruent'
    tmptmat(:,2:4) = tmat(:,1:3);
    tmptmat(:,1)   = 0;
    %[2 3] are incongruent
    tmptmat = tmptmat==2 | tmptmat==3;
    output  = sort(imat(tmptmat));
    output  = output(output>0);
  case 'previousleft'
  case 'previousright'
  otherwise
end
