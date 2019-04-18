clear all;


cd /project/3011085.03/analysis/granger

d=dir('*granger_pre.mat');
d1=dir('*prev1*');
d2=dir('*prev2*');
d3=dir('*prev3*');
d4=dir('*prev4*');
for k = 1:numel(d1)
  load(d(k).name);
  G(:,:,:,k) = granger.grangerspctrm;
  load(d1(k).name);
  G1(:,:,:,k) = granger.grangerspctrm;
  load(d2(k).name);
  G2(:,:,:,k) = granger.grangerspctrm;
  load(d3(k).name);
  G3(:,:,:,k) = granger.grangerspctrm;
  load(d4(k).name);
  G4(:,:,:,k) = granger.grangerspctrm;
end
dG13 = G1-G3;
dG42 = G4-G2;

%identify parcel combinations with zig-zags, due to spectral factorization
%convergence issues.
dX   = squeeze(mean(abs(diff(G,1,3)),3));
dX13 = squeeze(mean(abs(diff(dG13,1,3)),3));
dX42 = squeeze(mean(abs(diff(dG42,1,3)),3));

for k = 1:size(dX13,3)
  P = prctile(reshape(dX(:,:,k),[],1),[25 75]);
  outliers{k} = dX(:,:,k)>P(2)+4*diff(P);
  
  P = prctile(reshape(dX13(:,:,k),[],1),[25 75]);
  outliers13{k} = dX13(:,:,k)>P(2)+4*diff(P);
  P = prctile(reshape(dX42(:,:,k),[],1),[25 75]);
  outliers42{k} = dX42(:,:,k)>P(2)+4*diff(P);
end

for k = 1:19
  [x,y]=find(outliers{k});
  for m = 1:numel(x)
    G(x(m),y(m),:,k) = nan;
  end
  
  [x,y]=find(outliers13{k});
  for m = 1:numel(x)
    dG13(x(m),y(m),:,k) = nan;
    G1(x(m),y(m),:,k) = nan;
    G3(x(m),y(m),:,k) = nan;
    G4(x(m),y(m),:,k) = nan;
    G2(x(m),y(m),:,k) = nan;
  end
  [x,y]=find(outliers42{k});
  for m = 1:numel(x)
    dG42(x(m),y(m),:,k) = nan;
    G4(x(m),y(m),:,k) = nan;
    G2(x(m),y(m),:,k) = nan;
    G1(x(m),y(m),:,k) = nan;
    G3(x(m),y(m),:,k) = nan;
  end

end


G(~isfinite(G))= 0;
G1(~isfinite(G1))=0;
G2(~isfinite(G2))=0;
G3(~isfinite(G3))=0;
G4(~isfinite(G4))=0;

addpath ~/matlab/toolboxes/nmfv1_4
load vismot_parcels

P = P2*P;
n = size(P,1);

Gx = zeros(n,n,240,19);
G1x = zeros(n,n,240,19);
G2x = zeros(n,n,240,19);
G3x = zeros(n,n,240,19);
G4x = zeros(n,n,240,19);
for k = 1:19
  tmp  = G(:,:,1:end,k);
  tmp1 = G1(:,:,1:end,k);
  tmp2 = G2(:,:,1:end,k);
  tmp3 = G3(:,:,1:end,k);
  tmp4 = G4(:,:,1:end,k);
  
  for m = 1:size(tmp,3)
    Gx(:,:,m,k)  = P*tmp(:,:,m)*P';
    G1x(:,:,m,k) = P*tmp1(:,:,m)*P';
    G2x(:,:,m,k) = P*tmp2(:,:,m)*P';
    G3x(:,:,m,k) = P*tmp3(:,:,m)*P';
    G4x(:,:,m,k) = P*tmp4(:,:,m)*P';
       
  end
end

mGx = squeeze(nanmean(Gx,3));
for m = 1:size(mGx,3)
  qval = prctile(reshape(mGx(:,:,m),[],1),[25 75]);
  outliers_avg{m} = mGx(:,:,m)>1.5.*diff(qval)+qval(2);
end
% Gx2=Gx;
% for m= 1:size(Gx2,4)
%   tmp=Gx2(:,:,:,m);
%   for k = 1:240
%     tmptmp=tmp(:,:,k);
%     tmptmp(outliers_avg{m}|outliers_avg{m}')=nan;
%     tmp(:,:,k) = tmptmp;
%   end
%   Gx2(:,:,:,m)=tmp;
% end
% Gx2 = Gx;
% for m = 1:size(Gx2,4)
%   tmp = Gx2(:,:,:,m);
%   for k = 1:size(tmp,1)
%     for kk = (k+1):size(tmp,2)
%       tmptmp = squeeze(tmp([k kk],[k kk],:));
%       tmptmp = tmptmp./nansum(tmptmp);
%       tmp([k kk],[k kk],:) = tmptmp;
%     end
%   end
%   Gx2(:,:,:,m) = tmp;
% end
Gx2 = Gx;

dat = permute(Gx2,[3 4 1 2]);
datN = dat;
for k = 1:size(dat,2)
  %sdat(k) = nanstd(reshape(dat(:,k,:,:),[],1));
  %datN(:,k,:,:) = dat(:,k,:,:)./sdat(k);
  sdat(1,k,:,:) = nanmean(dat(:,k,:,:));
  datN(:,k,:,:) = dat(:,k,:,:)./sdat(1,k,:,:);
  %sdat2(k) = nanstd(reshape(datN(:,k,:,:),[],1));
  %datN(:,k,:,:) = datN(:,k,:,:)./sdat2(k);
end
datN(~isfinite(datN))=0;

% further reduce the size of the data
ulabel2 = {'R06';'R04';'R123';'R05';'R07';'R19';'R18';'R17';'L17';'L18';'L19';'L07';'L05';'L123';'L04';'L06'};
    
for k = 1:numel(ulabel2)
  indx = contains(ulabel,ulabel2{k});
  P2(k,indx) = 1./sum(indx);
end
dat2 = zeros(240,19,numel(ulabel2),numel(ulabel2));
for k = 1:240
  for m = 1:19
    tmp = shiftdim(dat(k,m,:,:));
    tmp2 = P2*tmp*P2';
    tmp2 = tmp2-diag(diag(tmp2));
    dat2(k,m,:,:) = tmp2;
  end
end
datN2 = dat2;
for k = 1:size(dat2,2)
  sdat(1,k,:,:) = nanmean(dat2(:,k,:,:));
  datN2(:,k,:,:) = dat2(:,k,:,:)./sdat(1,k,:,:);
end
datN2(~isfinite(datN2))=0;



option = [];
%option.eta = 200;
%option.beta = .8;
option.kernel = 'linear';
datN = reshape(datN, [240*19 n*n]);
indx = 1:n^2;
indx(1:(n+1):end) = [];
option.indx = indx;
[A,S]=sparsenmfnnls_wrapper(10,option);


nrand   = 50;
datadir = '/project/3011085.03/analysis/granger';

cd(datadir);
for k = 1:nrand
  % this will use the same randomseed for all # components
  option.randomseed = round(sum(100*clock));
  for ncomp = 5:30% (5:25)
    qsubfeval2('sparsenmfnnls_wrapper2',ncomp,option,'memreq',10*1024^3,'timreq',(30+ncomp)*60,'batchid',['nmf',num2str(k,'%03d'),'_ncomp',num2str(ncomp,'%03d')]);
  end
end

for m = 1:26
  ncomp = m + 4;
  
  d = dir(['nmf*ncomp',num2str(ncomp,'%03d'),'*output.mat']); nrand = numel(d);
  A = cell(1,nrand);
  S = cell(1,nrand);
  
  fprintf('processing ncomp=%d with %d randomizations\n', ncomp, nrand);
  for k = 1:nrand
    load(d(k).name);
    S{k} = argout{2};
    A{k} = argout{1};
    %ind  = getClusters(A{k}');
    ind = getClusters(S{k});
    if ~exist('cmat', 'var')
    %  cmat(size(A{1},1),size(A{1},1))=0;
      cmat(size(S{1},2),size(S{1},2))=0;    
    elseif exist('cmat', 'var') && k==1
      cmat(:)=0;
    end
    cmat = cmat+getRelationMatrix(ind);
    Res(m,k) = argout{5};
  end
  cmat = cmat./nrand;
  rho(m,1) = dispersionCoefficient(cmat);

  s = cat(1,S{:});
  %a{m} = cat(2,A{:});
  a{m} = cat(1,A{:})';
  
  allzero = sum(s==0,2)>size(s,2).*0.99;
  s(allzero,:) = [];
  a{m}(:,allzero) = [];
  
  C = s*s';
  C = C./sqrt(diag(C)*diag(C)');
  
  ft_hastoolbox('icasso',1);
  [P,Z,order] = hcluster(1-C, 'cl');
  R           = rindex(1-C, P);
  %[~,N]       = min(R(1:min(ncomp,numel(R))));
  [~,N]       = min(R(4:min(ncomp,numel(R)))); % JM changed 20190327 to avoid ending up in the 2-component case, which does not make sense really
  
  Rindex{m}   = R(1:min(ncomp,numel(R)));
  
  [iq,in,out] = clusterquality('mean',C,P(N,:));
  idx = zeros(N,1);
  for k = 1:N
    index        = find(P(N,:)==k);
    [~,idx(k,1)] = max(sum(C(index,index)));
    idx(k)       = index(idx(k));
  end
  
  [iq, ord] = sort(iq);
  Iq{m}     = iq;
  c{m}      = s(idx(ord),:);
  a{m}      = a{m}(:,idx(ord));
  a{m}      = reshape(a{m},[240 19 N]);
  idx = idx(ord);
end
save(fullfile(datadir,'groupresults_nmf'), 'c', 'a', 'Iq', 'Rindex', 'Res', 'rho');



% normalise the individual conditions subjectwise, based on the
% subject specific estimate of the data combined, but first NaN-out the
% autoconnections
for k = 1:size(Gx,1)
  Gx(k,k,:,:) = nan;
  G1x(k,k,:,:) = nan;
  G2x(k,k,:,:) = nan;
  G3x(k,k,:,:) = nan;
  G4x(k,k,:,:) = nan;
end

% M = nanmean(nanmean(nanmean(Gx,3),2),1);
% G1x = G1x./M;
% G2x = G2x./M;
% G3x = G3x./M;
% G4x = G4x./M;

label = cell(256,1);

for k = 1:256
  k1 = mod(k-1,16)+1;
  k2 = ceil(k/16);
  label{k} = sprintf('%s_%s',ulabel2{k1},ulabel2{k2});
end

% reorganize the data a bit
foi = 0:0.5:119.5;
data.freq = foi;
data.dimord = 'rpt_chan_freq';
data.label  = label;

data1 = data;
data2 = data;
data3 = data;
data4 = data;

data1.grangerspctrm = permute(reshape(G1x,[256 240 19]),[3 1 2]);
data2.grangerspctrm = permute(reshape(G2x,[256 240 19]),[3 1 2]);
data3.grangerspctrm = permute(reshape(G3x,[256 240 19]),[3 1 2]);
data4.grangerspctrm = permute(reshape(G4x,[256 240 19]),[3 1 2]);

data4f = data4;
data2f = data2;
data4f.grangerspctrm = permute(reshape(G4x(end:-1:1,end:-1:1,:),[256 240 19]),[3 1 2]);
data2f.grangerspctrm = permute(reshape(G2x(end:-1:1,end:-1:1,:),[256 240 19]),[3 1 2]);

datacmb1 = data1;
datacmb3 = data3;
datacmb1.grangerspctrm = (data1.grangerspctrm + data4f.grangerspctrm)./2;
datacmb3.grangerspctrm = (data3.grangerspctrm + data2f.grangerspctrm)./2;

n = size(G1x,4);
design = [ones(1,n) ones(1,n)*2;1:n 1:n];

cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 1000;
cfg.statistic = 'depsamplesT';
cfg.numrandomization = 1000;
cfg.correctm = 'cluster';
cfg.connectivity = speye(256);
cfg.clusteralpha = 0.01;
cfg.clusterthreshold = 'nonparametric_individual';
cfg.design = design;
cfg.method = 'montecarlo';
cfg.parameter = 'grangerspctrm';
cfg.neighbours = [];
%stat   = ft_freqstatistics(cfg, 
stat13 = ft_freqstatistics(cfg, data1,data3);
stat42f = ft_freqstatistics(cfg, data4f,data2f);
stat13cmb = ft_freqstatistics(cfg, datacmb1, datacmb3);


p13 = reshape(stat13.prob, [16 16 240]);
p42 = reshape(stat42.prob, [16 16 240]);
p   = reshape(stat.prob,   [16 16 240]);



