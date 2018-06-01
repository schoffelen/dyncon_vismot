function [trl,event] = trialfun_longtrials(cfg)

event = read_event(cfg.datafile);
value = [event.value]; %this can be shorter than event (because type occasionally is []
for k = 1:length(event)
  ok(k) = ~isempty(event(k).value);
end
event = event(ok);
value = value - bitand(value, 8192); %remove synchro trigger
value = value - bitand(value, 4096).*(4095./4096); %convert 4096 to 1
ok    = value>0;
value = value(ok);
event = event(ok);

trl = zeros(0,3);
for k = 1:length(event)-1
  if strcmp(event(k).type,'TRIGGER') && value(k)==8,
    cnt   = 0;
    kplus = 0;
    while cnt<4,
      kplus = kplus+1;
      if value(k+kplus)==32 && strcmp(event(k+kplus).type, 'TRIGGER'), cnt = cnt+1; end
    end
    smp    = event(k).sample;
    begsmp = event(k).sample;
    endsmp = event(k+kplus).sample;
    trlnew = [begsmp endsmp begsmp-smp+1]; %think about this
    trl    = [trl; trlnew];
  end
end
