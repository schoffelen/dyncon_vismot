function [vol] = prepareVol(mriname,name) 

load(mriname);
%mri = read_mri(mriname);
%cd('/raw/10/Anatomical/');
%cd([mriname,'_RAW']);
%d = dir;
%fname = d(3).name;
%
%%cd('/analyse/4/Project0030/pilot/');
%%fname = 'SubjectCMC.mri';
%mri   = read_fcdc_mri(fname);
%
%mriorig = mri;
%mriorig.transform1 = mriorig.transform;
%mriorig.transform = eye(4);
%clear mri;
%cfg     = [];
%cfg.method = 'interactive';
%mri     = volumerealign(cfg, mriorig);

cfg             = [];
cfg.coordinates = 'ctf';
cfg.template    = '/home/jan/matlab/spm2/templates/T1.mnc';
segment         = volumesegment(cfg, mri);
segment.gray    = flipdim(flipdim(flipdim(permute(segment.gray,[3 2 1]),3),2),1);
segment.white   = flipdim(flipdim(flipdim(permute(segment.white,[3 2 1]),3),2),1);
segment.csf     = flipdim(flipdim(flipdim(permute(segment.csf,[3 2 1]),3),2),1);
segment.dim     = size(segment.gray);


cfg = [];
vol = prepare_singleshell(cfg, segment);

fname = ['/analyse/4/Project0030/vol/',name,'vol'];
save(fname, 'vol');
