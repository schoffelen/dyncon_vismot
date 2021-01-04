function subject = vismot_subjinfo(subjectname)

if nargin<1
  subjectname = {'sub01';'sub02';'sub03';'sub04';...
                 'sub05';'sub06';'sub07';'sub08';...
                 'sub09';'sub10';'sub11';'sub12';...
                 'sub13';'sub14';'sub15';'sub16';...
                 'sub17';'sub18';'sub19'};
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
    praw = '/project/3011085.03/raw';
    pana = '/project/3011085.03/analysis';
end;
pwdir = pwd;

if ischar(subjectname)
  subjectname = {subjectname};
end

for k = 1:numel(subjectname)
  switch subjectname{k}
    case 'sub09'
      subject(k).name        = 'sub09';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'Seated1000/';
      subject(k).sessionname = '09-01-30@1414/';
      subject(k).runnames    = {'2/';'3/';'4/';'5/';'6/'};
      subject(k).emptyrun    = {'1/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'c,rfDC';
      subject(k).denoise     = 0;
      subject(k).congrROIfoi = [48 64; 40 52];
    case 'sub17'
      subject(k).name        = 'sub17';
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
    case 'sub12'
      subject(k).name        = 'sub12';
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
    case 'sub16'
      subject(k).name        = 'sub16';
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
    case 'sub08'
      subject(k).name        = 'sub08';
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
    case 'sub07'
      subject(k).name        = 'sub07';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-17@1042/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-17@1031/1/';
    case 'sub13'
      subject(k).name        = 'sub13';
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
    case 'sub15'
      subject(k).name        = 'sub15';
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
    case 'sub18'
      subject(k).name        = 'sub18';
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
    case 'sub06'
      subject(k).name        = 'sub06';
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
    case 'sub14'
      subject(k).name        = 'sub14';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-20@1310/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-20@1224/1/';
    case 'sub01'
      subject(k).name        = 'sub01';
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
    case 'sub19'
      subject(k).name        = 'sub19';
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
    case 'sub11'
      subject(k).name        = 'sub11';
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
    case 'sub05'
      subject(k).name        = 'sub05';
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
    case 'sub03'
      subject(k).name        = 'sub03';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-03-27@1613/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-03-27@1533/1/';
    case 'sub04'
      subject(k).name        = 'sub04';
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
    case 'sub02'
      subject(k).name        = 'sub02';
      subject(k).rawpath     = praw;
      subject(k).scanname    = 'S1k@6Ehmt/';
      subject(k).sessionname = '09-04-01@1036/';
      subject(k).runnames    = {'1/'};
      subject(k).emptyrun    = {'2/'};
      subject(k).pathname    = pana;
      subject(k).datafile    = 'hc,rfDC';
      subject(k).denoise     = 1;
      subject(k).cohname     = 'Std1k@6EEG/09-04-01@1028/1/';
      subject(k).ecgfile     = 'sub02ecg-aligned.mat';
      subject(k).ecgcomp     = [1 2 3];
      subject(k).congrROIfoi = [40 60; 64 80];
    case 'sub10'
      subject(k).name        = 'sub10';
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
