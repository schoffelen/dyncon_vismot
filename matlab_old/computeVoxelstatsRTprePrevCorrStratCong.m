function [stat13, stat42] = computeVoxelstatsRTpre(subject, frequency, smoothing)

%do sourcestatistics using precomputed fwhm and voxeldata 

%output:
%Descriptive statistics based on RT split
%stat13 stat42 divided according to response hand;

fieldtripdefs

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

cd([subject.pathname,'filter']);
if frequency ~= round(frequency)
  load([subject.name,'filt400_',num2str(round(10*smoothing),'%03d'),'_',num2str(round(10*frequency),'%04d')], 'fwhm');
else
  load([subject.name,'filt',num2str(round(10*smoothing),'%03d'),'_',num2str(frequency,'%03d')], 'fwhm');
end
cd([subject.pathname,'voxeldata']);
if frequency ~= round(frequency)
  load([subject.name,'voxeldata400_',num2str(round(10*smoothing),'%03d'),'_',num2str(round(10*frequency),'%04d')]);
else
  load([subject.name,'voxeldata',num2str(round(10*smoothing),'%03d'),'_',num2str(frequency,'%03d')]);
end
cd([subject.pathname,'grid']);
load([subject.name,'grid6mm.mat']);
eval('tmp = grid;');

%create symmetric inside
inside = isfinite(fwhm);% & isfinite(flipdim(fwhm,1));
inside = find(inside);
outside = setdiff(1:prod(tmp.dim), inside);
inside  = inside(:);
outside = outside(:);

%compute voxel specifice smoothing kernel
addpath /home/jan/projects/ccc/3D
tmp.fwhm    = fwhm;
tmp.inside  = inside;
tmp.outside = outside;
krn         = compute_kernel(tmp, 'truncate', 2e-5);
dim         = tmp.dim;

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.parameter = 'pow';
cfg2.statistic = 'indepsamplesT';
cfg2.ivar      = 1;

cfg3           = cfg2;
cfg3.statistic = 'yuent';

%congruency T-test
%activation condition 1 vs 3 response left
%activation condition 4 vs 2 response right
conditions = [1 3;4 2];
for mm = 1:2
  warning off
  tmprt = rtall;
  %exclude trials when previous trial was incorrect
  sel   = selectTrials(subject, 'previouscorrect');
  
  indx1 = intersect(btrl(conditions(mm,1)):etrl(conditions(mm,1)), sel(:)');
  indx2 = intersect(btrl(conditions(mm,2)):etrl(conditions(mm,2)), sel(:)');
  indx  = [indx1 indx2];
  
  selPC  = selectTrials(subject, 'previouscongruent');
  selPiC = selectTrials(subject, 'previousincongruent');
  selPN  = selectTrials(subject, 'previousneutral');

  %median split
  rt    = tmprt(indx);
  rt    = rt(:)';
  mrt   = median(rt); %determine median of only the candidate trials
  irt   = iqr(rt);
  thrlo = mrt - 1.5*irt;
  thrhi = mrt + 1.5*irt;
  tmprt(rtall<thrlo) = nan;
  tmprt(rtall>thrhi) = nan;
  
  dummy = prctile(rt,[33.3333333 66.6666667]);
    
  selslow = find(tmprt>=dummy(2)); %relative to all trials
  selfast = find(tmprt< dummy(1));

  slowPC  = intersect(intersect(selslow(:), selPC(:)) , indx(:));;
  slowPiC = intersect(intersect(selslow(:), selPiC(:)), indx(:));;
  slowPN  = intersect(intersect(selslow(:), selPN(:)) , indx(:));;
  fastPC  = intersect(intersect(selfast(:), selPC(:)) , indx(:));;
  fastPiC = intersect(intersect(selfast(:), selPiC(:)), indx(:));;
  fastPN  = intersect(intersect(selfast(:), selPN(:)) , indx(:));;

  minPC  = min(numel(slowPC), numel(fastPC));
  minPiC = min(numel(slowPiC),numel(fastPiC));
  minPN  = min(numel(slowPN), numel(fastPN));
  
  dum = randperm(minPC) ;slowPC  = sort(slowPC(dum(1:minPC)));
  dum = randperm(minPiC);slowPiC = sort(slowPiC(dum(1:minPiC)));
  dum = randperm(minPN) ;slowPN  = sort(slowPN(dum(1:minPN)));
  dum = randperm(minPC) ;fastPC  = sort(fastPC(dum(1:minPC)));
  dum = randperm(minPiC);fastPiC = sort(fastPiC(dum(1:minPiC)));
  dum = randperm(minPN) ;fastPN  = sort(fastPN(dum(1:minPN)));

  indx1 = sort([fastPC;fastPiC;fastPN]); n1 = numel(indx1);
  indx2 = sort([slowPC;slowPiC;slowPN]); n2 = numel(indx2);
  indx  = [indx1(:)' indx2(:)'];

  tmptrl = findcfg(sd.cfg,'trl');
  if size(tmptrl,1) ~= size(sd.pow,1),
    sd.cfg.trl = zeros(size(sd.pow,1),3);
  end
  tmpsd = selectdata(sd, 'rpt', indx);
  
  tmpsd.pow = standardise(log10(tmpsd.pow), 1);
  cfg2.design      = [ones(1,n1) ones(1,n2)*2]; 
  cfg3.design      = cfg2.design;

  tmpstat         = sourcestatistics(cfg2, tmpsd);
  tmpstat.inside  = inside;
  tmpstat.outside = outside;
  tmpstat.pos     = tmp.pos;
  dum             = zeros(size(tmpstat.pos,1),1);
  dum(inside)     = tmpstat.stat;
  tmpstat.stat    = dum;
  tmpstat         = rmfield(tmpstat,'mask');
  tmpstat         = rmfield(tmpstat,'prob');
  %tmpstat         = rmfield(tmpstat,'dim');
  tmpstat.stat2   = smooth_vol(tmpstat.stat,krn,dim,inside);
  tmpstat.cfg     = [];
  tmpstaty        = sourcestatistics(cfg3, tmpsd);
  dum(inside)     = tmpstaty.stat;
  tmpstat.staty   = dum;
  tmpstat.stat2y  = smooth_vol(tmpstat.staty,krn,dim,inside); 
  
  if mm==1,
    stat13 = tmpstat;
  elseif mm==2,
    stat42 = tmpstat;
  end
end
