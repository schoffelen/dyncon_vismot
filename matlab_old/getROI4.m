function [roi1,roi2,inside,roivol,foi] = getROI4(subject, saveflag, foilim)

if nargin<2,
  saveflag = 1;
end

datdir = '/analyse/4/Project0030/source/4D';
cd(datdir);

load([subject.name,'stat4D']);

dim = pos2dim3d(stat13.pos);
dum = zeros(dim);

sel1 = nearest(stat13.freq, foilim(1,1)):nearest(stat13.freq, foilim(1,2));
sel2 = nearest(stat13.freq, foilim(2,1)):nearest(stat13.freq, foilim(2,2));
dum  = reshape(mean(stat13.stat2(:,sel1),2),dim);
dum1 = tfce(dum,[],[],[],median(dum(stat13.inside)));
dum  = reshape(mean(stat42.stat2(:,sel2),2),dim);
if strcmp(subject.name, 'T001'),
  dum2 = tfce(dum);
else
  dum2 = tfce(dum,[],[],[],median(dum(stat13.inside)));
end
mask = zeros(dim);
mask(1:15,1:15,11:18) = 1;
mask(17:31,1:15,11:18) = 2;
mask(stat13.outside) = 0;

head    = zeros(dim);
head(:) = 1:prod(dim);
head(stat13.outside) = nan;
left    = head;
left(mask~=1) = nan;
right   = head;
right(mask~=2) = nan;
ix{1} = left(find(isfinite(left)));
ix{2} = right(find(isfinite(right)));
roivol = zeros(dim);
for k = 1:2
  
  if k==1,
    tmpvol = double(isfinite(left)).*dum2;
  else
    tmpvol = double(isfinite(right)).*dum1;
  end
  roivol = tmpvol + roivol;
end
%tmpvol = double(isfinite(left)).*-1.*reshape(mean(stat.stat3y(:,foilim),2),dim);
%roivol = tmpvol;
%tmpvol = double(isfinite(right)).*-1.*reshape(mean(stat.stat3y(:,foilim),2),dim);
%roivol = tmpvol + roivol;


ok = [0 0];

blob1 = roivol.*double(mask==1);
blob2 = roivol.*double(mask==2);

[srt1,indx1] = sort(blob1(:),'descend');
[srt2,indx2] = sort(blob2(:),'descend');

while sum(ok)~=2,

  roi1 = indx1(1:25);
  roi2 = indx2(1:25);

  dum = zeros(dim);
  dum(roi1)=1;
  dum = bwlabeln(dum);
  if max(dum(:))>1,
    for j=1:max(dum(:))
      n(j) = sum(dum(:)==j);
    end
    [mn,in] = max(n);
    sel     = find(dum>0 & dum~=in);
    indx1(ismember(indx1,sel)) = [];
  else
    ok(1)=1;
  end
  
  dum = zeros(dim);
  dum(roi2)=1;
  dum = bwlabeln(dum);
  if max(dum(:))>1,
    for j=1:max(dum(:))
      n(j) = sum(dum(:)==j);
    end
    [mn,in] = max(n);
    sel     = find(dum>0 & dum~=in);
    indx2(ismember(indx2,sel)) = [];
  else
    ok(2)=1;
  end

end

inside = stat13.inside;
foi    = stat13.freq([sel1([1 end]);sel2([1 end])]);
%foi = [8 10 12];
if saveflag,
  cd('/analyse/4/Project0030/roi');
  save([subject.name,'roimaskCongruency'],'roi1','roi2','inside','roivol','foi');
end
