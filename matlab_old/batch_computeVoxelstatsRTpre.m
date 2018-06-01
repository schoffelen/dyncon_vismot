subjinfo;

%savepath = '/analyse/4/Project0030/voxelstats';
savepath = '/analyse/1/Project0002/tmpProject0030/voxelstats';
smoothing = 3.2;
freqs     = 22.4;
freqno    = 1;

for k = 1:numel(selno)
  subject = SUBJ(selno(k));
  for m = 1:numel(freqno)
    frequency    = freqs(freqno(m));
    [stat13, stat42] = computeVoxelstatsRTpre(subject, frequency, smoothing);
    %[stat13, stat42] = computeVoxelstatsRTprePreviouscorrect(subject, frequency, smoothing);
    %[stat13, stat42] = computeVoxelstatsRTprePrevCorrStratCong(subject, frequency, smoothing);
    
    warning off;
    stat13 = struct2single(stat13);
    stat42 = struct2single(stat42);

    cd(savepath);
    %save([subject.name,'voxelstats400_',num2str(round(10*smoothing),'%03d'), ...
    %     '_',num2str(round(10*frequency),'%04d')], 'stat13','stat42');
    save([subject.name,'voxelstats625_',num2str(round(10*smoothing),'%03d'), ...
         '_',num2str(round(10*frequency),'%03d')], 'stat13','stat42');
    clear stat13 stat42
  end
end
