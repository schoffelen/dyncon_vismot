foi = [16:2:28];
for k = 1:7
  fstr   = num2str(foi(k),'%03d');
  suffix = ['stat',fstr,'pre2roicong_038'];
  [statRi,statRc,statLi,statLc,roi1,roi2,statri,statrc,statli,statlc] = collectVoxelstatsCOH2('/analyse/4/Project0030/source/',suffix,1,0);
  figure;
  subplot(2,2,2);volplotJM(reshape(statRc.stat,dim),'montage');title(['roi1r ',fstr]);drawnow;
  subplot(2,2,3);volplotJM(reshape(statLi.stat,dim),'montage');title(['roi1l ',fstr]);drawnow;
  subplot(2,2,4);volplotJM(reshape(statRi.stat,dim),'montage');title(['roi2r ',fstr]);drawnow;
  subplot(2,2,1);volplotJM(reshape(statLc.stat,dim),'montage');title(['roi2l ',fstr]);drawnow;
end

foi = [16:2:28];
roi1l = zeros(15*7,33480);
roi2l = zeros(15*7,33480);
roi1r = zeros(15*7,33480);
roi2r = zeros(15*7,33480);
for k = 1:7
  fstr   = num2str(foi(k),'%03d');
  suffix = ['stat',fstr,'pre2roicong_038'];
  [statRi,statRc,statLi,statLc,roi1,roi2,statri,statrc,statli,statlc] = collectVoxelstatsCOH2('/analyse/4/Project0030/source/',suffix,1,0);
  roi1l((k-1)*15+[1:15],:) = statli.dcoh(1:15,:);
  roi2l((k-1)*15+[1:15],:) = statlc.dcoh(1:15,:);
  roi1r((k-1)*15+[1:15],:) = statrc.dcoh(1:15,:);
  roi2r((k-1)*15+[1:15],:) = statri.dcoh(1:15,:);
end
x1l = corr(roi1l(:,inside)');
x2l = corr(roi2l(:,inside)');
x1r = corr(roi1r(:,inside)');
x2r = corr(roi2r(:,inside)');

%sort frequencies per subject
ix = [];
for k = 1:15
  ix = [ix k:15:105];
end
x1l = x1l(ix,ix);
x2l = x2l(ix,ix);
x1r = x1r(ix,ix);
x2r = x2r(ix,ix);

for k = 1:15
  x1l((k-1)*7+[1:7],(k-1)*7+[1:7]) = 0;
  x2l((k-1)*7+[1:7],(k-1)*7+[1:7]) = 0;
  x1r((k-1)*7+[1:7],(k-1)*7+[1:7]) = 0;
  x2r((k-1)*7+[1:7],(k-1)*7+[1:7]) = 0;
end

tmpfoi=4;
figure;
for k = 1:15
  volplotJM(reshape(roi1r((k-1)*15+tmpfoi,:),dim),'montage');pause;clf;
end

ix = reshape(1:105, [15 7])';
for k = 1:15
  ixmax1(:,k) = max(max(roi1l(ix(:,k),:),[],2));
  ixmin1(:,k) = min(min(roi1l(ix(:,k),:),[],2));
  ixmax2(:,k) = max(max(roi2r(ix(:,k),:),[],2));
  ixmin2(:,k) = min(min(roi2r(ix(:,k),:),[],2));
end
crange1 = max(abs(ixmin1),ixmax1);
crange2 = max(abs(ixmin2),ixmax2);

for k = 1:15
  figure;
  for m = 1:7
    ix = mod(m-1,3);
    iy = floor((m-1)/3);
    subplot('position',[0.35*ix 0.7-0.35*iy 0.3 0.3]);
    volplotJM(reshape(roi1l((m-1)*15+k,:),dim)+head,'montage');title([num2str(k),' ',num2str(m)]); caxis([-1 1].*crange1(k));
  end
  dum = zeros(dim);
  dum(roi1{k}) = 1;
  dum(statri.outside) = -1;
  subplot('position',[0.35 0 0.3 0.3]);
  volplotJM(dum,'montage');
  print(gcf,'-dpng',['subject',num2str(k,'%03d'),'roi1l']);
end

