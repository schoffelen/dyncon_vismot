subjinfo;

sel = find(~ismember({SUBJ(:).name},{'GAR12' 'BKA01' 'KBI24'}));
for k = 1:numel(sel)
  dat(:,:,k) = SUBJ(sel(k)).rs.nmin1trim(:,1:4);
end
ix = size(dat,1);
iy = size(dat,2);
iz = size(dat,3);

%subject
tmp      = zeros(1,1,iz);
tmp(:)   = 1:iz;
tmp      = repmat(tmp, [ix iy 1]);
group{1} = tmp(:);

%congruency trial n
tmp           = zeros(ix,iy);
tmp([1 4], :) = 1;
tmp([2 3], :) = 2;
tmp           = repmat(tmp, [1 1 iz]);
group{2}      = tmp(:);

%response side trial n
tmp           = zeros(ix,iy);
tmp([1 3], :) = 1;
tmp([2 4], :) = 2;
tmp           = repmat(tmp, [1 1 iz]);
group{3}      = tmp(:);

%congruency trial n-1
tmp           = zeros(ix,iy);
tmp(:,[1 4])  = 1;
tmp(:,[2 3])  = 2;
tmp           = repmat(tmp, [1 1 iz]);
group{4}      = tmp(:);

%response side trial n-1
tmp           = zeros(ix,iy);
tmp(:,[1 3])  = 1;
tmp(:,[2 4])  = 2;
tmp           = repmat(tmp, [1 1 iz]);
%group{4}      = tmp(:);
group{5}      = tmp(:);

groupnames = {'subject' 'congruencyN' 'responsesideN' 'congruencyNmin1' 'responseNmin1'};
%groupnames = {'subject' 'congruencyN' 'responsesideN' 'congruencyNmin1'};
%groupnames = {'subject' 'congruencyN' 'responsesideN' 'responseNmin1'};
[p,t,stats,terms] = anovan(dat(:), group, 'varnames', groupnames, 'random', 1);
[p2,t2,stats2,terms2] = anovan(dat(:), group, 'varnames', groupnames, 'model', 'interaction', 'random', 1);
[p3,t3,stats3,terms3] = anovan(dat(:), group, 'varnames', groupnames, 'model', 3, 'random', 1);

[c,m,h,g] = multcompare(stats3, 'dimension', [1 2], 'display', 'off');
mc12.c    = c;
mc12.m    = m;
mc12.h    = h;
mc12.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [1 3], 'display', 'off');
mc13.c    = c;
mc13.m    = m;
mc13.h    = h;
mc13.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [1 4], 'display', 'off');
mc14.c    = c;
mc14.m    = m;
mc14.h    = h;
mc14.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [1 5], 'display', 'off');
mc15.c    = c;
mc15.m    = m;
mc15.h    = h;
mc15.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 3], 'display', 'off');
mc23.c    = c;
mc23.m    = m;
mc23.h    = h;
mc23.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 4], 'display', 'off');
mc24.c    = c;
mc24.m    = m;
mc24.h    = h;
mc24.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 5], 'display', 'off');
mc25.c    = c;
mc25.m    = m;
mc25.h    = h;
mc25.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [3 4], 'display', 'off');
mc34.c    = c;
mc34.m    = m;
mc34.h    = h;
mc34.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [3 5], 'display', 'off');
mc35.c    = c;
mc35.m    = m;
mc35.h    = h;
mc35.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [4 5], 'display', 'off');
mc45.c    = c;
mc45.m    = m;
mc45.h    = h;
mc45.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 3 4], 'display', 'off');
mc234.c    = c;
mc234.m    = m;
mc234.h    = h;
mc234.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 3 5], 'display', 'off');
mc235.c    = c;
mc235.m    = m;
mc235.h    = h;
mc235.g    = g;
[c,m,h,g] = multcompare(stats3, 'dimension', [2 4 5], 'display', 'off');
mc245.c    = c;
mc245.m    = m;
mc245.h    = h;
mc245.g    = g;

x = reshape(mc24.m(:,1)./1.01725,[2 2]);% x = x(:,[1 3 2]); %put neutral in middle
figure;bar(x); title('CongNvsCongNmin1');
axis([0.5 2.5 500 700]);

x = reshape(mc23.m(:,1)./1.01725,[2 2]);
figure;bar(x'); title('CongNvsRespside');
axis([0.5 2.5 500 700]);

x = reshape(mc234.m([1 3 5 7 2 4 6 8],1)./1.01725,[4 2]);
figure;bar(x'); title('CongNvsRespsideNvsCongNmin1');
axis([0.5 2.5 500 700]);

x = reshape(mc235.m([1 5 3 7 2 6 4 8],1)./1.01725,[4 2]);
figure;bar(x'); title('CongNvsRespsideNvsRespsideNmin1');
axis([0.5 2.5 500 700]);

x = reshape(mc245.m([1 5 3 7 2 6 4 8],1)./1.01725,[4 2]);
figure;bar(x'); title('CongNvsCongNmin1vsRespsideNmin1');
axis([0.5 2.5 500 700]);
