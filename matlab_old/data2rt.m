function [rt] = data2rt(subject, flag)

%compute rt for a given subject
%flag refers to data_aligned (1) or data(0)

if nargin<2, flag = 0; end

fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1},subject.datafile];
hdr     = read_header(fname);
fsample = hdr.Fs;

if length(subject.runnames)>1,
  for k = 1:length(subject.runnames)
    fname   = [subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{k},subject.datafile];
    event{k} = ft_read_event(fname);
  end
else
  event{1} = ft_read_event(fname);
end

cd(subject.pathname);
cd('data');
if flag,
  load([subject.name,'data_aligned']);
else
  load([subject.name,'data']);
end

for k = 1:5
  warning off;
  if k==1,
    data = struct2double(data1);
  elseif k==2,
    data = struct2double(data2);
  elseif k==3,
    data = struct2double(data3);
  elseif k==4,
    data = struct2double(data4);
  elseif k==5,
    data = struct2double(data5);
  end
  warning on;
  
  if k<6,
    trl = findcfg(data.cfg, 'trl');
    if length(subject.runnames)>1,
      %make cell-array out of trl corresponding with the runs
      tmp  = diff([trl(:,1);0]);
      btrl = [0;find(tmp(1:end-1)<0)]+1;
      etrl = [find(tmp<0)];
      for kk = 1:numel(btrl)
        newtrl{kk,1} = trl(btrl(kk):etrl(kk),:);
      end
      trl = newtrl;
    else
      trl   = {trl};
    end

    rt{k} = [];
    
    for kk = 1:length(event)
      if k~=5,
        cfg = [];
        cfg.eventtype   = 'RESPONSE';
        cfg.eventvalue  = 2.^(3*mod(k+1,2)+5);
        cfg.searchrange = 'afterzero';
        cfg.output      = 'samplefromoffset';
        tmprt = ft_recodeevent(cfg, event{kk}, trl{kk});
        tmprt = tmprt(:);
      else
        tmprt = zeros(size(trl{kk},1),1)+nan;
      end

      tmp = subject.trigsmp;
      tmp(subject.runnr(:,1)~=kk) = nan;
      for mm = 1:size(trl{kk},1)
        [ix,iy]    = find(tmp == trl{kk}(mm,1)-trl{kk}(mm,3)+1);
        tmprt(mm,2) = subject.trlid(ix,iy); 
      end
      
      rt{k} = [rt{k}; tmprt]; %in number of samples
    end
  else
    rt{k} = nan+zeros(numel(data.trial),2); 
  end
end
