subjinfo;

savepath = '/analyse/1/Project0002/tmpProject0030/filter';
%freqs    = [1:1:100];
smoothing = 3.2;
freqs  = 22.4;
freqno = 1;
selno  = find(~ismember({SUBJ(:).name}, {'BKA01' 'GAR12' 'KBI24'}));

for k = 1:numel(selno)
%for k = 12:numel(selno)
  subject = SUBJ(selno(k))
  for m = 1:numel(freqno)
    frequency    = freqs(freqno(m));
    [filt, fwhm] = computeFiltersPre(subject, frequency, smoothing);
    cd(savepath);
    %save([subject.name,'filt',num2str(round(10*smoothing),'%03d'), ...
    %     '_',num2str(frequency,'%03d')], 'filt', 'fwhm');
    save([subject.name,'filt625_',num2str(round(10*smoothing),'%03d'), ...
         '_',num2str(frequency,'%03d')], 'filt', 'fwhm');
    clear filt fwhm
  end
end
