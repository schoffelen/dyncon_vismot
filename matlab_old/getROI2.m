function [roi1,roi2,inside] = getROI2(subject)

datdir = '/analyse/4/Project0030/source/4Dprepst';
cd(datdir);

load([subject.name,'stat4Dprepstbroadband']);

dim = pos2dim3d(stat.pos);
foi = [3:5 7:12 18:23];
fac = 1;

%subject T001 has no nice gamma response, take beta instead
%if strcmp(subject.name,'T001') || strcmp(subject.name,'TMR04')
%  foi = 7:12;
%  fac = -1;
%end
%tmp = fac.*reshape(mean(stat.stat2(:,foi),2),dim);
%get rois for posterior blobs of induced gamma

tmp = reshape(mean(1./stat.fwhm(:,foi),2), dim);
tmp = tfce(tmp);

mask = zeros(dim);
mask(1:15,1:15,11:18) = 1;
mask(17:31,1:15,11:18) = 2;
mask(stat.outside) = 0;

ok = [0 0];

tmp(~isfinite(tmp))=0;
blob1 = tmp.*double(mask==1);
blob2 = tmp.*double(mask==2);

[srt1,indx1] = sort(blob1(:),'descend');
[srt2,indx2] = sort(blob2(:),'descend');

while sum(ok)~=2,
  roi1 = indx1(1:40);
  roi2 = indx2(1:40);

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

cd('/analyse/4/Project0030/roi');
save([subject.name,'roimaskFWHMabg'],'roi1','roi2','inside');

