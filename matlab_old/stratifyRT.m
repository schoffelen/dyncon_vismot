function [output, input, binaxis] = stratifyRT(subject,  data)

if nargin<2,
  data = [];
end

if ~iscell(data) && (isempty(data) || data==0),
  cd(subject.pathname);
  cd('rt');
  load([subject.name,'rt']);
else
  for k = 1:4
    tmp = findcfg(data{k}.cfg, 'trl');
    rt{k} = tmp(:,4);
  end
end

cfg               = [];
cfg.equalbinavg   = 'no';
cfg.numbin        = 8;
input{1,1}        = rt{1}';
input{1,2}        = rt{2}';
input{1,3}        = rt{3}';
input{1,4}        = rt{4}';
[output, binaxis] = stratify(cfg,input);
