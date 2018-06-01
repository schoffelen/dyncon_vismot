cd('/raw/11/Project0030/mss15/');

%first get orientation and positions of coils
%in dewar space, using the COH-measurement
%these data will serve as a constraint for the
%fitting procedure in the next step
cd('Std1k@6EEG/09-04-03@1052/1');
hdr = read_header('c,rfDC');
ubs = hdr.orig.user_block_data;
coh = ubs{15}; %lucky guess, first one is usually 15, second one 17


cd('/raw/11/Project0030/mss15/');
cd('Test@Movin/09-04-03@1104/1/');
event = read_event('hc,rfDC');
type  = {event.type}';
sel   = find(strmatch('TRIGGER', type));
val   = [event(sel).value]';
smp   = [event(sel).sample]';

block = find(val==256);
for k = 1:(length(block)-1)
  tmpval = val(block(k)+1:block(k+1)-1);
  tmpsmp = smp(block(k)+1:block(k+1)-1);
  sel    = find(tmpval==4224);
  trl{k} = [tmpsmp(sel)-500 tmpsmp(sel)+499 ones(length(sel),1)-501];
  trl{k} = trl{k}(2:end-1,:);
end

for k = 1:length(trl)
  cfg          = [];
  cfg.datafile = 'hc,rfDC';
  cfg.channel  = {'MEG' 'UACurrent'};
  cfg.trl      = trl{k};
  cfg.padding  = 5;
  cfg.hpfilter = 'yes';
  cfg.hpfreq   = 10;
  data         = preprocessing(cfg);
  
  cfg = [];
  cfg.coils.pos = coh.pnt;
  cfg.coils.ori = coh.ori;
  cfg.twindow   = 'trial';
  hmdata{k}     = headmotiontracking(cfg, data);
end

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
    plot3(pnt{m}(1:2:248,1,k),pnt{m}(1:2:248,2,k),pnt{m}(1:2:248,3,k),'.','markersize',4);
  end
end



