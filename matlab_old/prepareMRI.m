function [mri] = prepareMRI(mriname, shapename) 

cd('/raw/10/Anatomical/');
cd([mriname,'_RAW']);
d = dir;
fname = d(3).name;

mri   = read_mri(fname);
shape = read_headshape(shapename);

cfg   = [];
mri   = mycoreg(cfg,mri,shape);
