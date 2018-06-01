subjinfo;

savepath = '/analyse/4/Project0030/voxelstats';
freqs    = [1:1:100];
smoothing = 3.75;

for k = 1:numel(selno)
  subject = SUBJ(selno(k));
  for m = 1:numel(freqno)
    frequency    = freqs(freqno(m));
    [stat13, stat42] = computeVoxelstatsPCpre(subject, frequency, smoothing);
    
    warning off;
    stat13 = struct2single(stat13);
    stat42 = struct2single(stat42);

    cd(savepath);
    save([subject.name,'voxelstats',num2str(round(10*smoothing),'%03d'), ...
         '_',num2str(frequency,'%03d')], 'stat13','stat42');
    clear stat13 stat42
  end
end
