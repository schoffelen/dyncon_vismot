function [roi1,roi2,inside,roivol,foi,stat] = getROI3(subject, saveflag, foilim)

if nargin<2,
  saveflag = 1;
end

datdir = '/analyse/4/Project0030/source/4Dprepst';
cd(datdir);

load([subject.name,'stat4Dprepstbroadband']);

if nargin<3,
  foilim = 1:numel(stat.freq);
  %foilim = 1:12;
end

dim = pos2dim3d(stat.pos);
dum = zeros(dim);
for k = foilim
  fprintf('tfce %d/%d\n', k, numel(stat.freq));
  dum = reshape(stat.stat2y(:,k), dim);
  dum = tfce(dum);
  stat.stat3y(:,k) = dum(:);
end

mask = zeros(dim);
mask(1:15,1:15,11:18) = 1;
mask(17:31,1:15,11:18) = 2;
mask(stat.outside) = 0;

head    = zeros(dim);
head(:) = 1:prod(dim);
head(stat.outside) = nan;
left    = head;
left(mask~=1) = nan;
right   = head;
right(mask~=2) = nan;
ix{1} = left(find(isfinite(left)));
ix{2} = right(find(isfinite(right)));
roivol = zeros(dim);
for k = 1:2
  
  tmpmax = max(stat.stat3y(ix{k},foilim),[],1);
  tmpmin = min(stat.stat3y(ix{k},foilim),[],1);
  [m, iy] = max([tmpmax; abs(tmpmin)], [], 1);
  foi(k) = find(m==max(m));
  if iy(foi(k))==1,
    fac = 1;
  else
    fac = -1;
  end
  if k==1,
    tmpvol = double(isfinite(left)).*fac.*reshape(mean(stat.stat3y(:,foilim(foi(k))),2),dim);
  else
    tmpvol = double(isfinite(right)).*fac.*reshape(mean(stat.stat3y(:,foilim(foi(k))),2),dim);
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

inside = stat.inside;
foi    = stat.freq(foilim(foi));
%foi = [8 10 12];
if saveflag,
  cd('/analyse/4/Project0030/roi');
  save([subject.name,'roimaskTFCE'],'roi1','roi2','inside','roivol','foi');
end
