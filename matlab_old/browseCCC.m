foi = [16:2:28];
for k = 1:7
  fstr   = num2str(foi(k),'%03d');
  suffix = ['stat',fstr,'pre2roicong_038'];
  [stat,statl,statr,roi1,roi2,stat1l,stat2l,stat1r,stat2r,stat1,stat2] = collectVoxelstatsCOH('/analyse/4/Project0030/source/',suffix,1,0);
  %figure;volplotJM(reshape(stat1r.stat,dim),'montage');title(['roi1r ',fstr]);drawnow;
  %figure;volplotJM(reshape(stat2l.stat,dim),'montage');title(['roi2l ',fstr]);drawnow;
  %figure;volplotJM(reshape(stat2r.stat,dim),'montage');title(['roi2r ',fstr]);drawnow;
  %figure;volplotJM(reshape(stat1l.stat,dim),'montage');title(['roi1l ',fstr]);drawnow;
  figure;
  subplot(2,2,2);volplotJM(reshape(stat1r.stat,dim),'montage');title(['roi1r ',fstr]);drawnow;
  subplot(2,2,3);volplotJM(reshape(stat2l.stat,dim),'montage');title(['roi2l ',fstr]);drawnow;
  subplot(2,2,4);volplotJM(reshape(stat2r.stat,dim),'montage');title(['roi2r ',fstr]);drawnow;
  subplot(2,2,1);volplotJM(reshape(stat1l.stat,dim),'montage');title(['roi1l ',fstr]);drawnow;
end

foi = [16:2:28];
roi1l = zeros(15*7,33480);
roi2l = zeros(15*7,33480);
roi1r = zeros(15*7,33480);
roi2r = zeros(15*7,33480);
for k = 1:7
  fstr   = num2str(foi(k),'%03d');
  suffix = ['stat',fstr,'pre2roicong_038'];
  [stat,statl,statr,roi1,roi2,stat1l,stat2l,stat1r,stat2r,stat1,stat2] = collectVoxelstatsCOH('/analyse/4/Project0030/source/',suffix,1,0);
  roi1l((k-1)*15+[1:15],:) = statl.dcoh1s(1:15,:);
  roi2l((k-1)*15+[1:15],:) = statl.dcoh2s(1:15,:);
  roi1r((k-1)*15+[1:15],:) = statr.dcoh1s(1:15,:);
  roi2r((k-1)*15+[1:15],:) = statr.dcoh2s(1:15,:);
end
x1l = corr(roi1l(:,inside)');
x2l = corr(roi2l(:,inside)');
x1r = corr(roi1r(:,inside)');
x2r = corr(roi2r(:,inside)');

%fort frequencies per subject
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

tmpfoi=6;
figure;
for k = 1:15
  volplotJM(reshape(roi1r((k-1)*15+tmpfoi,:),dim),'montage');pause;clf;
end
