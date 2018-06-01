function [data1,data2,data3,data4,data5] = doRealignment(subject)

cd(subject.pathname);
cd('data');
load([subject.name,'data']);
cd('../dewar2head_avg');
load([subject.name,'dewar2head_avg']);
cd('../vol');
vol = read_vol([subject.name,'vol.mat']);

for k = 1:5
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
  warning off;
  data            = struct2double(data);
  data.grad       = convert_units(data.grad, 'cm');
  cfg             = [];
  %cfg.template{1} = transform_sens(M, convert_units(data.grad, 'cm'));
  %cfg.vol      = vol;
  cfg.template{1} = data.grad;
  cfg.feedback      = 'no';
  cfg.vol         = transform_vol(inv(M), vol);
  cfg.inwardshift = 1;
  data            = megrealign(cfg, data);
  data            = struct2single(data);
  if k==1,
    data1 = data;
  elseif k==2,
    data2 = data;
  elseif k==3,
    data3 = data;
  elseif k==4,
    data4 = data;
  elseif k==5,
    data5 = data;
  end
end
