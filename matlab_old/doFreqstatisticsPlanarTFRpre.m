function [stat13, stat42] = doFreqstatisticsPlanarTFRpre(subject)

%do sourcestatistics using precomputed fwhm and voxeldata 

%output:
%Descriptive statistics based on RT split
%stat13 stat42 divided according to response hand;

fieldtripdefs

cd([subject.pathname,'freq']);
load([subject.name,'tfrpreHanning']);

cfg2           = [];
cfg2.method    = 'montecarlo';
cfg2.numrandomization = 0;
cfg2.parameter = 'powspctrm';
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

  %exclude trials when previous trial was incorrect
  sel     = selectTrials(subject, 'previouscorrect', allfreq(conditions(mm,:)));
  tmpfreq = selectdata(allfreq{conditions(mm,:)}, 'param', 'powspctrm'); 
  rtall   = [];
  dummy   = findcfg(allfreq{conditions(mm,1)}.cfg, 'trl');
  rtall   = [rtall;dummy(:,4)];
  dummy   = findcfg(allfreq{conditions(mm,2)}.cfg, 'trl');
  rtall   = [rtall;dummy(:,4)];

  %%median split
  %rt    = freq.rt(indx);
  %rt    = rt(:)';
  %mrt   = median(rt);
  %rt(rt<mrt)  = 1; %fast
  %rt(rt>=mrt) = 2; %slow
  %tmpsd = selectdata(sd, 'rpt', indx);
 
  tmptrl = findcfg(tmpfreq.cfg,'trl');
  if size(tmptrl,1) ~= size(tmpfreq.powspctrm,1),
    tmpfreq.cfg.trl = zeros(size(tmpfreq.powspctrm,1),3);
  end
  tmpsd = selectdata(tmpfreq, 'rpt', sel);
  tmprt = rtall(sel);

  %exclude outliers >1.5 iqr from median
  tmprt = tmprt(:)';
  mrt   = median(tmprt);
  irt   = iqr(tmprt);
  thrlo = mrt - 1.5*irt;
  thrhi = mrt + 1.5*irt;
  tmprt(tmprt<thrlo)=nan;
  tmprt(tmprt>thrhi)=nan;
  sel   = find(isfinite(tmprt));
  tmprt = tmprt(sel);
  tmpfreq = selectdata(tmpfreq, 'rpt', sel);
  [srt,srtix] = sort(tmprt);
  %nx    = ceil(numel(tmprt)/3);
  nx    = ceil(numel(tmprt)/4);
  indx  = srtix([1:nx (numel(srtix)-nx+1):numel(srtix)]);
  tmprt = tmprt(indx);
  tmpfreq.cfg.trl = zeros(size(tmpfreq.powspctrm,1),3);
  tmpfreq = selectdata(tmpfreq, 'rpt', indx);
  rt    = [ones(1,nx) ones(1,nx)*2];

  tmpfreq.powspctrm = nanstandardise(log10(tmpfreq.powspctrm), 1);
  cfg2.design      = rt; 
  cfg3.design      = cfg2.design;

  tmpstat         = freqstatistics(cfg2, tmpfreq);
  dum     = tmpstat.stat;
  dum(isnan(dum)) = 0;
  dum             = spm_t2z(dum, size(cfg2.design,2)-2);
  tmpstat.stat    = dum;
  tmpstaty        = freqstatistics(cfg3, tmpfreq);
  dum     = tmpstaty.stat;
  dum(isnan(dum)) = 0;
  dum             = spm_t2z(dum, size(cfg2.design,2)-2);
  tmpstat.staty   = dum;
  
  if mm==1,
    stat13 = tmpstat;
  elseif mm==2,
    stat42 = tmpstat;
  end
end
