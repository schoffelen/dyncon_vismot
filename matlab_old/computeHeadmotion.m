function [data] = computeHeadmotion(data, hmtlist, subject)

if nargin==1,
  hmtlist = 'MEG';
end

if ~isempty(subject.cohname),
  %first get orientation and positions of coils
  %in dewar space, using the COH-measurement
  %these data will serve as a constraint for the
  %fitting procedure in the next step
  cd([subject.rawpath,subject.name,'/',subject.cohname]);
  hdr = read_header('c,rfDC');
  ubs = hdr.orig.user_block_data;
  coh1 = ubs{15}; %lucky guess, first one is usually 15, second one 17
  if strcmp(subject.name, 'BKA01'),
    fprintf('swapping coil order\n');
    coh1.pnt = coh1.pnt([2 1 3 5 4],:);
    coh1.ori = coh1.ori([2 1 3 5 4],:);
  end
  coh2 = ubs{12}; %lucky guess, digitized head points
else
  %there is no first guess, analyze some data to get there
  cd([subject.rawpath,subject.name,'/',subject.scanname,subject.sessionname,subject.runnames{1}]);
  hdr  = read_header('hc,rfDC');
  ubs = hdr.orig.user_block_data;
  coh2 = ubs{12};
  
  cfg          = [];
  cfg.datafile = 'hc,rfDC';
  cfg.channel  = [hmtlist; {'UACurrent'}];
  cfg.trl      = round([hdr.Fs*10+1 hdr.Fs*20 0]);
  cfg.padding  = 15;
  cfg.hpfilter = 'yes';
  cfg.hpfreq   = 10;
  tmpdata      = preprocessing(cfg);
  
  cfg           = [];
  cfg.coils.dewar.pos  = warp_apply([eye(3) [0;0;-0.14];[0 0 0 1]],coh2.pnt(6:10,:));
  cfg.coils.head.pos   = coh2.pnt;
  cfg.coils.head.label = coh2.label;
  cfg.twindow   = 'trial';
  cfg.channel   = tmpdata.label(1:end-1);
  cfg.resample  = 'no';
  [tmphmdata,tmpori] = headmotiontracking(cfg, tmpdata);
  dat           = tmphmdata.trial{1};
  M             = headcoordinates(dat(1:3)',dat(4:6)',dat(7:9)');
  coh1.pnt      = warp_apply(inv(M), coh2.pnt(6:10,:));
  %coh1.ori      = [];
  coh1.ori      = tmpori;
end

nrun = length(subject.runnames);

alltrl       = findcfg(data.cfg, 'trl');;
if nrun>1,
  drun = [find(diff([inf;alltrl(:,1)])<0);size(alltrl,1)+1];
  for k = 1:length(drun)-1
    trl{k} = alltrl(drun(k):drun(k+1)-1,:);
    ntrl(1,k) = size(trl{k},1);
  end
else
  trl{1} = alltrl;
  ntrl   = size(trl{1},1);
end

for k = 1:nrun
  cd([subject.rawpath,subject.name,'/',subject.scanname,...
      subject.sessionname,subject.runnames{k}]);

  cfg          = [];
  cfg.datafile = 'hc,rfDC';
  cfg.channel  = [hmtlist; {'UACurrent'}];
  cfg.trl      = trl{k};
  %cfg.trl = cfg.trl(1:5,:);
  cfg.padding  = 5;
  cfg.hpfilter = 'yes';
  cfg.hpfreq   = 10;
  tmpdata      = preprocessing(cfg);
  
  cfg           = [];
  cfg.coils.dewar.pos  = coh1.pnt;
  if isfield(coh1, 'ori'), cfg.coils.dewar.ori  = coh1.ori; end
  cfg.coils.head.pos   = coh2.pnt;
  cfg.coils.head.label = coh2.label;
  cfg.twindow   = 'trial';
  cfg.channel   = tmpdata.label(1:end-1);
  cfg.resample  = 'no';
  hmdata{k}     = headmotiontracking(cfg, tmpdata);
end

%resample to the input data
ntrl = [0 cumsum(ntrl)];
for k = 1:length(hmdata)
  cfg         = [];
  cfg.detrend = 'no';
  cfg.time    = data.time(ntrl(k)+1:ntrl(k+1));
  hmdata{k}   = resampledata(cfg, hmdata{k});
end

if length(hmdata)>1,
  tmpdata = appenddata([], hmdata{:});
  hmdata    = {tmpdata};
end

%FIXME assume 1 run only
data = appenddata([], data, hmdata{1});

doplot = 0;
if doplot,

grad = hmdata{1}.grad;
for m = 1:length(hmdata)
  fprintf('analysing block %d/%d\n',m,length(hmdata));
  pnt{m} = zeros([size(grad.pnt) length(hmdata{m}.trial)]);
  for k = 1:length(hmdata{m}.trial)
    M = rigidbodyJM(hmdata{m}.trial{k}(1:6));
    tmpgrad = transform_sens(M, grad);
    pnt{m}(:,:,k) = tmpgrad.pnt;
  end
end

for m = 1:length(pnt)
  figure;hold on;
  for k = 1:size(pnt{m},3);
    plot3(pnt{m}(1:248,1,k),pnt{m}(1:248,2,k),pnt{m}(1:248,3,k),'.','markersize',4);
  end
end

mpnt     = mean(pnt{1}(1:248,:,:),3);
npnt     = size(mpnt,1);
gradorig = hdr.grad;
pntorig  = gradorig.pnt;

T = ([pntorig(1:248,:) ones(npnt,1)]')/([mpnt ones(npnt,1)]');
T(1,1:3) = T(1,1:3)./norm(T(1,1:3));
T(2,1:3) = T(2,1:3)./norm(T(2,1:3));
T(3,1:3) = T(3,1:3)./norm(T(3,1:3));
T(4,:)   = [0 0 0 1];

end

