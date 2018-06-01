subjinfo;

savepath = '/analyse/1/Project0002/tmpProject0030/voxeldata';
smoothing = 3.2;
freqs  = 22.4;
freqno = 1;
selno  = find(~ismember({SUBJ(:).name}, {'BKA01' 'GAR12' 'KBI24'}));

for k = 1:numel(selno)
%for k = 2:numel(selno)
  subject = SUBJ(selno(k));
  for m = 1:numel(freqno)
    frequency    = freqs(freqno(m));
    [sd,btrl,etrl,rtall] = computeVoxeldata(subject, frequency, smoothing);
    cd(savepath);
    %save([subject.name,'voxeldata',num2str(round(10*smoothing),'%03d'), ...
    %     '_',num2str(frequency,'%03d')], 'sd', 'btrl','etrl','rtall');
    save([subject.name,'voxeldata625_',num2str(round(10*smoothing),'%03d'), ...
         '_',num2str(round(10*frequency),'%03d')], 'sd', 'btrl','etrl','rtall');
    clear sd btrl etrl rtall
  end
end
