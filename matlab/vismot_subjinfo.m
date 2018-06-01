function subject = vismot_subjinfo(subjectname)

if nargin<1
  subjectname = {'AHE08';'AMN20';'BKA01';'CBE22';...
                 'CDE04';'DBD12';'GAR12';'JGN27';...
                 'JOE22';'LCE09';'MAI27';'MHH14';...
                 'MME25';'PCL19';'RBE13';'SZQ02';...
                 'T001';'TMR04';'VIA12'};
end

if strcmp(computer, 'PCWIN')
    praw = 'U:/';
    pana = 'Y:/';
elseif strcmp(computer, 'MACI'),
    praw = '/Volumes/11/11/';
    pana = '/Volumes/ANALYSE/4/';
    %come up with something
else
    %praw = '/raw/11/';
    %pana = '/analyse/4/';
    praw = '/home/language/jansch/projects/visuomotor/data/raw';
    pana = '/home/language/jansch/projects/visuomotor/data/analyse';
end;
pwdir = pwd;

if ischar(subjectname)
  subjectname = {subjectname};
end

for k = 1:numel(subjectname)
  switch subjectname{k}
    case 'JOE22'
      subject(k).name        = 'JOE22';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'Seated1000/';
      subject(k).sessionname = '09-01-30@1414/';
      subject(k).runnames    = {'2/';'3/';'4/';'5/';'6/'};
      subject(k).emptyrun    = {'1/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'c,rfDC';
      subject(k).denoise     = 0;
      subject(k).congrROIfoi = [48 64; 40 52];
    case 'T001'
      subject(k).name        = 'T001';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'Seated1000/';
      subject(k).sessionname = '09-02-05@1014/';
      subject(k).runnames    = {'3/';'4/'};
      subject(k).emptyrun    = {'1/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'c,rfDC';
      subject(k).mriname     = 'GGA04';
      subject(k).denoise     = 0;
      subject(k).congrROIfoi = [36 64; 84 96];
    case 'MHH14'
      subject(k).name        = 'MHH14';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'SVeNoeeg10/';
      subject(k).sessionname = '09-02-19@1225/';
      subject(k).runnames    = {'1/';'2/'};
      subject(k).emptyrun    = {'3/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'c,rfDC';
      subject(k).denoise     = 0;
      subject(k).congrROIfoi = [20 20; 44 56];
    case 'KBI24'
      subject(k).name        = 'KBI24';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'JMSFinal/';
      subject(k).sessionname = '09-03-03@1447/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
    case 'SZQ02'
      subject(k).name        = 'SZQ02';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'Seated1Khm/';
      subject(k).sessionname = '09-03-06@1409/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      %subject(k).badchannels = {'-A107';'-A121';'-A40';'-A146';'-A248'};
      subject(k).denoise     = 1;
      subject(k).congrROIfoi = [20 20; 60 80];
    case 'JGN27'
      subject(k).name        = 'JGN27';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1kchmt/';
      subject(k).sessionname = '09-03-10@1105/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      %subject(k).badchannels = {'-A107';'-A40';'-A248'};
      subject(k).denoise     = 1;
      subject(k).congrROIfoi = [64 88; 72 92];
    case 'GAR12'
      subject(k).name        = 'GAR12';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-17@1042/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-17@1031/1/';
    case 'MME25'
      subject(k).name        = 'MME25';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-18@1634/';
      subject(k).runnames    = {'1/' '2/'};
      subject(k).emptyrun    = {'3/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-18@1621/1/';
      subject(k).congrROIfoi = [60 80; 60 80];
    case 'RBE13'
      subject(k).name        = 'RBE13';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-19@1206/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-19@1201/1/';
      subject(k).congrROIfoi = [60 80; 60 80];
    case 'TMR04'
      subject(k).name        = 'TMR04';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-19@1500/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-19@1450/1/';
      subject(k).congrROIfoi = [52 60; 44 52];
    case 'DBD12'
      subject(k).name        = 'DBD12';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-20@1026/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-20@0940/1/';
      subject(k).congrROIfoi = [40 60; 40 60];
    case 'PCL19'
      subject(k).name        = 'PCL19';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-20@1310/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-20@1224/1/';
    case 'AHE08'
      subject(k).name        = 'AHE08';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-20@1605/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-20@1524/1/';
      subject(k).congrROIfoi = [60 80; 60 80];
    case 'VIA12'
      subject(k).name        = 'VIA12';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-23@1201/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-23@1112/1/';
      subject(k).congrROIfoi = [80 92; 68 88];
    case 'MAI27'
      subject(k).name        = 'MAI27';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-25@1338/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-25@1334/1/';
      subject(k).congrROIfoi = [48 60; 60 80];
    case 'CDE04'
      subject(k).name        = 'CDE04';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-25@1603/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-25@1510/1/';
      subject(k).congrROIfoi = [36 52; 36 56];
    case 'BKA01'
      subject(k).name        = 'BKA01';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-27@1613/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-27@1533/1/';
    case 'CBE22'
      subject(k).name        = 'CBE22';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-31@1043/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-31@1034/1/';
      subject(k).congrROIfoi = [40 60; 36 52];
    case 'AMN20'
      subject(k).name        = 'AMN20';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-04-01@1036/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-04-01@1028/1/';
      subject(k).ecgfile     = 'AMN20ecg-aligned.mat';
      subject(k).ecgcomp     = [1 2 3];
      subject(k).congrROIfoi = [40 60; 64 80];
    case 'LCE09'
      subject(k).name        = 'LCE09';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-04-01@1311/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-04-01@1304/1/';
      subject(k).congrROIfoi = [36 56; 68 80];
    otherwise
      error('unknown subject requested');
  end
end

for k = 1:numel(subject)
  tmp = subject(k);
  if ~exist(fullfile(tmp.pathname,'behaviour',[tmp.name,'behaviour.mat']), 'file')
    [rt, trl, correct, trigsmp, rs, runnr, trlid, startfix] = analyzeRT2(tmp);
    save(fullfile(tmp.pathname,'behaviour',[tmp.name,'behaviour']), 'rt', 'trl', 'correct', 'trigsmp', 'rs', 'trlid', 'runnr', 'startfix');
  else
    load(fullfile(tmp.pathname,'behaviour',[tmp.name,'behaviour']));
  end
  subject(k).rt      = rt;
  subject(k).trl     = trl;
  subject(k).correct = correct;
  subject(k).trigsmp = trigsmp;
  subject(k).runnr   = runnr;
  subject(k).trlid   = trlid;
  subject(k).startfix = startfix;
  subject(k).rs      = rs;
   
  clear rt trl correct rs trigsmp;
end
