function [data1,data2,data3,data4,data5] = cleanECGaligned(subject,flag)

if nargin==1,
  %compute decomposition flag
  flag = 1;
end

%load the data from which the spatial components will be removed
cd([subject.pathname,'data/']);

load([subject.name,'data_aligned']);
warning off
data = struct2double(appenddata([], data1,data2,data3,data4,data5));
warning on

%check whether data.cfg contains a resampletrl
if isempty(findcfg(data.cfg, 'resampletrl')),
  data = fixResample(data);
else
  %do nothing
end

if flag,
  compecg = clean_ecg2(data);
  cd([subject.pathname,'ecg/']);
  save([subject.name,'ecg-aligned'], 'compecg');
  
  cfg        = [];
  cfg.layout = '4D248.lay';
  componentbrowser(cfg, compecg);
  x = str2num(input('which components to remove? ','s'));
  subject.ecgcomp = x;
else
  keyboard
  %FIXME implement this
  %load the file with the spatial components to be removed
  cd([subject.pathname,'ecg/']);
  load(subject.ecgfile);
end

%get intersection of channels present in all conditions
alllabel    = data.label;
cfg         = [];
cfg.channel = alllabel;
clear data;
warning off;
data1 = struct2double(preprocessing(cfg, data1)); 
data2 = struct2double(preprocessing(cfg, data2)); 
data3 = struct2double(preprocessing(cfg, data3)); 
data4 = struct2double(preprocessing(cfg, data4)); 
data5 = struct2double(preprocessing(cfg, data5)); 
warning on;

cfg           = [];
cfg.component = subject.ecgcomp;
warning off;
data1 = struct2single(rejectcomponent(cfg,compecg,data1));
data2 = struct2single(rejectcomponent(cfg,compecg,data2));
data3 = struct2single(rejectcomponent(cfg,compecg,data3));
data4 = struct2single(rejectcomponent(cfg,compecg,data4));
data5 = struct2single(rejectcomponent(cfg,compecg,data5));

