function [s] = collect4D(subjname, savesuffix, loadsuffix, foi);

if nargin<2,
  savesuffix = '';
end
if nargin<3,
  loadsuffix = 'cohPreCongClean';
end
if nargin<4,
  foi = [8:2:36 40:4:100];
end

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);

for k = 1:numel(foi)
  k
  fname = [subjname,'stat',num2str(foi(k),'%03d'),loadsuffix];
  load(fname);
  for kk = 1:4
    if kk==1,
      stat = statlc;
    elseif kk==2,
      stat = statli;
    elseif kk==3,
      stat = statrc;
    elseif kk==4,
      stat = statri;
    end
    stat.freq = foi(k);
    stat.coh1 = trimmean(stat.coh1,0.2,2);
    stat.coh2 = trimmean(stat.coh2,0.2,2);
    stat.coh3 = trimmean(stat.coh3,0.2,2);
    stat.coh4 = trimmean(stat.coh4,0.2,2);

    if k==1,
      s{kk,1} = stat;
      s{kk,1}.coh1(:,2:numel(foi)) = 0;
      s{kk,1}.coh2(:,2:numel(foi)) = 0;
      s{kk,1}.coh3(:,2:numel(foi)) = 0;
      s{kk,1}.coh4(:,2:numel(foi)) = 0;
      s{kk,1}.freq(2:numel(foi)) = 0;
      s{kk,1}.fwhm(:,2:numel(foi)) = 0;
    else
      s{kk,1}.coh1(:,k) = stat.coh1;
      s{kk,1}.coh2(:,k) = stat.coh2;
      s{kk,1}.coh3(:,k) = stat.coh3;
      s{kk,1}.coh4(:,k) = stat.coh4;
      s{kk,1}.freq(k)   = foi(k);
      s{kk,1}.fwhm(:,k) = stat.fwhm(:);
    end
  end
end
%s{1}.dcoh = s{1}.coh1-s{1}.coh2;
%s{2}.dcoh = s{2}.coh1-s{2}.coh2;
%s{3}.dcoh = s{3}.coh1-s{3}.coh2;
%s{4}.dcoh = s{4}.coh1-s{4}.coh2;

if ~isempty(savesuffix)
  save([subjname,savesuffix], 's');
end
