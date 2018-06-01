function [s] = collect4Dpow3(subjname, savesuffix, loadsuffix, foi);

if nargin<2,
  savesuffix = '';
end
if nargin<3,
  loadsuffix = 'cohPreCongClean';
end
if nargin<4,
  foi = 8:2:26;
end

datadir = '/data1/synchro1/Projects/JanMathijs/Project0030tmp/source';
cd(datadir);

for k = 1:numel(foi)
  k
  fname = [subjname,'stat',num2str(foi(k),'%03d'),loadsuffix];
  load(fname);
  for kk = 1:2
    if kk==1,
      stat = statl;
    elseif kk==2,
      stat = statr;
    end
    stat.freq = foi(k);

    if k==1,
      s{kk,1} = stat;
      s{kk,1}.stat1(:,2:numel(foi)) = 0;
      s{kk,1}.stat2(:,2:numel(foi)) = 0;
      s{kk,1}.stat3(:,2:numel(foi)) = 0;
      s{kk,1}.stat4(:,2:numel(foi)) = 0;
      s{kk,1}.freq(2:numel(foi)) = 0;
      s{kk,1}.fwhm(:,2:numel(foi)) = 0;
    else
      s{kk,1}.stat1(:,k) = stat.stat1;
      s{kk,1}.stat2(:,k) = stat.stat2;
      s{kk,1}.stat3(:,k) = stat.stat3;
      s{kk,1}.stat4(:,k) = stat.stat4;
      s{kk,1}.freq(k)   = foi(k);
      s{kk,1}.fwhm(:,k) = stat.fwhm(:);
    end
  end
end

if ~isempty(savesuffix)
  save([subjname,savesuffix], 's');
end
