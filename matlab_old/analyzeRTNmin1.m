subjinfo;

sel = find(~ismember({SUBJ(:).name},{'GAR12' 'BKA01' 'KBI24'}));
for k = 1:numel(sel)
  dat(:,:,k) = SUBJ(sel(k)).rs.nmin1trim;
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
tmp(:,[1 4])  = -1;
tmp(:,[2 3])  = 0;
tmp(:,5)      = 1;
tmp           = repmat(tmp, [1 1 iz]);
group{4}      = tmp(:);

%response side trial n-1
tmp           = zeros(ix,iy);
tmp(:,[1 3])  = -1;
tmp(:,[2 4])  = 0;
tmp(:,5)      = 1;
tmp           = repmat(tmp, [1 1 iz]);
%group{4}      = tmp(:);
group{5}      = tmp(:);

groupnames = {'subject' 'congruencyN' 'responsesideN' 'congruencyNmin1' 'responseNmin1'};
%groupnames = {'subject' 'congruencyN' 'responsesideN' 'congruencyNmin1'};
%groupnames = {'subject' 'congruencyN' 'responsesideN' 'responseNmin1'};
[p,t,stats,terms] = anovan(dat(:), group, 'varnames', groupnames, 'random', 1);
[p2,t2,stats2,terms2] = anovan(dat(:), group, 'varnames', groupnames, 'model', 'interaction', 'random', 1);
[p3,t3,stats3,terms3] = anovan(dat(:), group, 'varnames', groupnames, 'model', 3, 'random', 1);

[c,m,h,g] = multcompare(stats2, 'dimension', [1 2], 'display', 'off');
mc12.c    = c;
mc12.m    = m;
mc12.h    = h;
mc12.g    = g;
[c,m,h,g] = multcompare(stats2, 'dimension', [1 3], 'display', 'off');
mc13.c    = c;
mc13.m    = m;
mc13.h    = h;
mc13.g    = g;
[c,m,h,g] = multcompare(stats2, 'dimension', [1 4], 'display', 'off');
mc14.c    = c;
mc14.m    = m;
mc14.h    = h;
mc14.g    = g;
[c,m,h,g] = multcompare(stats2, 'dimension', [2 3], 'display', 'off');
mc23.c    = c;
mc23.m    = m;
mc23.h    = h;
mc23.g    = g;
[c,m,h,g] = multcompare(stats2, 'dimension', [2 4], 'display', 'off');
mc24.c    = c;
mc24.m    = m;
mc24.h    = h;
mc24.g    = g;
[c,m,h,g] = multcompare(stats2, 'dimension', [3 4], 'display', 'off');
mc34.c    = c;
mc34.m    = m;
mc34.h    = h;
mc34.g    = g;



x = reshape(mc24.m(:,1)./1.01725,[2 3]); x = x(:,[1 3 2]); %put neutral in middle
figure;bar(x); title('CongNvsCongNmin1');
axis([0.3 2.7 500 700]);

x = reshape(mc23.m(:,1)./1.01725,[2 2]);
figure;bar(x'); title('CongNvsRespside');
axis([0.5 2.5 500 700]);


