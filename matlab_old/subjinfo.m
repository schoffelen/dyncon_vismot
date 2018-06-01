if strmatch(computer, 'PCWIN')
    praw = 'U:/';
    pana = 'Y:/';
elseif strmatch(computer, 'MACI'),
    praw = '/Volumes/11/11/';
    pana = '/Volumes/ANALYSE/4/';
    %come up with something
else
    praw = '/raw/11/';
    pana = '/analyse/4/';
end;
pwdir = pwd;

SUBJ(1).name        = 'JOE22';
SUBJ(1).rawpath     = [praw,'Project0030/'];
SUBJ(1).scanname    = 'Seated1000/';
SUBJ(1).sessionname = '09-01-30@1414/';
SUBJ(1).runnames    = {'2/';'3/';'4/';'5/';'6/'};
SUBJ(1).emptyrun    = {'1/'};
SUBJ(1).pathname    = [pana,'Project0030/'];
SUBJ(1).datafile    = 'c,rfDC';
SUBJ(1).denoise     = 0;
SUBJ(1).congrROIfoi = [48 64; 40 52];

SUBJ(2).name        = 'T001';
SUBJ(2).rawpath     = [praw,'Project0030/'];
SUBJ(2).scanname    = 'Seated1000/';
SUBJ(2).sessionname = '09-02-05@1014/';
SUBJ(2).runnames    = {'3/';'4/'};
SUBJ(2).emptyrun    = {'1/'};
SUBJ(2).pathname    = [pana,'Project0030/'];
SUBJ(2).datafile    = 'c,rfDC';
SUBJ(2).mriname     = 'GGA04';
SUBJ(2).denoise     = 0;
SUBJ(2).congrROIfoi = [36 64; 84 96];

SUBJ(3).name        = 'MHH14';
SUBJ(3).rawpath     = [praw,'Project0030/'];
SUBJ(3).scanname    = 'SVeNoeeg10/';
SUBJ(3).sessionname = '09-02-19@1225/';
SUBJ(3).runnames    = {'1/';'2/'};
SUBJ(3).emptyrun    = {'3/'};
SUBJ(3).pathname    = [pana,'Project0030/'];
SUBJ(3).datafile    = 'c,rfDC';
SUBJ(3).denoise     = 0;
SUBJ(3).congrROIfoi = [20 20; 44 56];

SUBJ(4).name        = 'KBI24';
SUBJ(4).rawpath     = [praw,'Project0030/'];
SUBJ(4).scanname    = 'JMSFinal/';
SUBJ(4).sessionname = '09-03-03@1447/';
SUBJ(4).runnames    = {'1/'};
SUBJ(4).emptyrun    = {'2/'};
SUBJ(4).pathname    = [pana,'Project0030/'];
SUBJ(4).datafile    = 'hc,rfDC';
SUBJ(4).denoise     = 1;

SUBJ(5).name        = 'SZQ02';
SUBJ(5).rawpath     = [praw,'Project0030/'];
SUBJ(5).scanname    = 'Seated1Khm/';
SUBJ(5).sessionname = '09-03-06@1409/';
SUBJ(5).runnames    = {'1/'};
SUBJ(5).emptyrun    = {'2/'};
SUBJ(5).pathname    = [pana,'Project0030/'];
SUBJ(5).datafile    = 'hc,rfDC';
%SUBJ(5).badchannels = {'-A107';'-A121';'-A40';'-A146';'-A248'};
SUBJ(5).denoise     = 1;
SUBJ(5).congrROIfoi = [20 20; 60 80];

SUBJ(6).name        = 'JGN27';
SUBJ(6).rawpath     = [praw,'Project0030/'];
SUBJ(6).scanname    = 'S1kchmt/';
SUBJ(6).sessionname = '09-03-10@1105/';
SUBJ(6).runnames    = {'1/'};
SUBJ(6).emptyrun    = {'2/'};
SUBJ(6).pathname    = [pana,'Project0030/'];
SUBJ(6).datafile    = 'hc,rfDC';
%SUBJ(6).badchannels = {'-A107';'-A40';'-A248'};
SUBJ(6).denoise     = 1;
SUBJ(6).congrROIfoi = [64 88; 72 92];

SUBJ(7).name        = 'GAR12';
SUBJ(7).rawpath     = [praw,'Project0030/'];
SUBJ(7).scanname    = 'S1k@6Ehmt/';
SUBJ(7).sessionname = '09-03-17@1042/';
SUBJ(7).runnames    = {'1/'};
SUBJ(7).emptyrun    = {'2/'};
SUBJ(7).pathname    = [pana,'Project0030/'];
SUBJ(7).datafile    = 'hc,rfDC';
SUBJ(7).denoise     = 1;
SUBJ(7).cohname     = 'Std1k@6EEG/09-03-17@1031/1/';

SUBJ(8).name        = 'MME25';
SUBJ(8).rawpath     = [praw,'Project0030/'];
SUBJ(8).scanname    = 'S1k@6Ehmt/';
SUBJ(8).sessionname = '09-03-18@1634/';
SUBJ(8).runnames    = {'1/' '2/'};
SUBJ(8).emptyrun    = {'3/'};
SUBJ(8).pathname    = [pana,'Project0030/'];
SUBJ(8).datafile    = 'hc,rfDC';
SUBJ(8).denoise     = 1;
SUBJ(8).cohname     = 'Std1k@6EEG/09-03-18@1621/1/';
SUBJ(8).congrROIfoi = [60 80; 60 80];

SUBJ(9).name        = 'RBE13';
SUBJ(9).rawpath     = [praw,'Project0030/'];
SUBJ(9).scanname    = 'S1k@6Ehmt/';
SUBJ(9).sessionname = '09-03-19@1206/';
SUBJ(9).runnames    = {'1/'};
SUBJ(9).emptyrun    = {'2/'};
SUBJ(9).pathname    = [pana,'Project0030/'];
SUBJ(9).datafile    = 'hc,rfDC';
SUBJ(9).denoise     = 1;
SUBJ(9).cohname     = 'Std1k@6EEG/09-03-19@1201/1/';
SUBJ(9).congrROIfoi = [60 80; 60 80];

SUBJ(10).name        = 'TMR04';
SUBJ(10).rawpath     = [praw,'Project0030/'];
SUBJ(10).scanname    = 'S1k@6Ehmt/';
SUBJ(10).sessionname = '09-03-19@1500/';
SUBJ(10).runnames    = {'1/'};
SUBJ(10).emptyrun    = {'2/'};
SUBJ(10).pathname    = [pana,'Project0030/'];
SUBJ(10).datafile    = 'hc,rfDC';
SUBJ(10).denoise     = 1;
SUBJ(10).cohname     = 'Std1k@6EEG/09-03-19@1450/1/';
SUBJ(10).congrROIfoi = [52 60; 44 52];

SUBJ(11).name        = 'DBD12';
SUBJ(11).rawpath     = [praw,'Project0030/'];
SUBJ(11).scanname    = 'S1k@6Ehmt/';
SUBJ(11).sessionname = '09-03-20@1026/';
SUBJ(11).runnames    = {'1/'};
SUBJ(11).emptyrun    = {'2/'};
SUBJ(11).pathname    = [pana,'Project0030/'];
SUBJ(11).datafile    = 'hc,rfDC';
SUBJ(11).denoise     = 1;
SUBJ(11).cohname     = 'Std1k@6EEG/09-03-20@0940/1/';
SUBJ(11).congrROIfoi = [40 60; 40 60];

SUBJ(12).name        = 'PCL19';
SUBJ(12).rawpath     = [praw,'Project0030/'];
SUBJ(12).scanname    = 'S1k@6Ehmt/';
SUBJ(12).sessionname = '09-03-20@1310/';
SUBJ(12).runnames    = {'1/'};
SUBJ(12).emptyrun    = {'2/'};
SUBJ(12).pathname    = [pana,'Project0030/'];
SUBJ(12).datafile    = 'hc,rfDC';
SUBJ(12).denoise     = 1;
SUBJ(12).cohname     = 'Std1k@6EEG/09-03-20@1224/1/';

SUBJ(13).name        = 'AHE08';
SUBJ(13).rawpath     = [praw,'Project0030/'];
SUBJ(13).scanname    = 'S1k@6Ehmt/';
SUBJ(13).sessionname = '09-03-20@1605/';
SUBJ(13).runnames    = {'1/'};
SUBJ(13).emptyrun    = {'2/'};
SUBJ(13).pathname    = [pana,'Project0030/'];
SUBJ(13).datafile    = 'hc,rfDC';
SUBJ(13).denoise     = 1;
SUBJ(13).cohname     = 'Std1k@6EEG/09-03-20@1524/1/';
SUBJ(13).congrROIfoi = [60 80; 60 80];

SUBJ(14).name        = 'VIA12';
SUBJ(14).rawpath     = [praw,'Project0030/'];
SUBJ(14).scanname    = 'S1k@6Ehmt/';
SUBJ(14).sessionname = '09-03-23@1201/';
SUBJ(14).runnames    = {'1/'};
SUBJ(14).emptyrun    = {'2/'};
SUBJ(14).pathname    = [pana,'Project0030/'];
SUBJ(14).datafile    = 'hc,rfDC';
SUBJ(14).denoise     = 1;
SUBJ(14).cohname     = 'Std1k@6EEG/09-03-23@1112/1/';
SUBJ(14).congrROIfoi = [80 92; 68 88];

SUBJ(15).name        = 'MAI27';
SUBJ(15).rawpath     = [praw,'Project0030/'];
SUBJ(15).scanname    = 'S1k@6Ehmt/';
SUBJ(15).sessionname = '09-03-25@1338/';
SUBJ(15).runnames    = {'1/'};
SUBJ(15).emptyrun    = {'2/'};
SUBJ(15).pathname    = [pana,'Project0030/'];
SUBJ(15).datafile    = 'hc,rfDC';
SUBJ(15).denoise     = 1;
SUBJ(15).cohname     = 'Std1k@6EEG/09-03-25@1334/1/';
SUBJ(15).congrROIfoi = [48 60; 60 80];

SUBJ(16).name        = 'CDE04';
SUBJ(16).rawpath     = [praw,'Project0030/'];
SUBJ(16).scanname    = 'S1k@6Ehmt/';
SUBJ(16).sessionname = '09-03-25@1603/';
SUBJ(16).runnames    = {'1/'};
SUBJ(16).emptyrun    = {'2/'};
SUBJ(16).pathname    = [pana,'Project0030/'];
SUBJ(16).datafile    = 'hc,rfDC';
SUBJ(16).denoise     = 1;
SUBJ(16).cohname     = 'Std1k@6EEG/09-03-25@1510/1/';
SUBJ(16).congrROIfoi = [36 52; 36 56];

SUBJ(17).name        = 'BKA01';
SUBJ(17).rawpath     = [praw,'Project0030/'];
SUBJ(17).scanname    = 'S1k@6Ehmt/';
SUBJ(17).sessionname = '09-03-27@1613/';
SUBJ(17).runnames    = {'1/'};
SUBJ(17).emptyrun    = {'2/'};
SUBJ(17).pathname    = [pana,'Project0030/'];
SUBJ(17).datafile    = 'hc,rfDC';
SUBJ(17).denoise     = 1;
SUBJ(17).cohname     = 'Std1k@6EEG/09-03-27@1533/1/';

SUBJ(18).name        = 'CBE22';
SUBJ(18).rawpath     = [praw,'Project0030/'];
SUBJ(18).scanname    = 'S1k@6Ehmt/';
SUBJ(18).sessionname = '09-03-31@1043/';
SUBJ(18).runnames    = {'1/'};
SUBJ(18).emptyrun    = {'2/'};
SUBJ(18).pathname    = [pana,'Project0030/'];
SUBJ(18).datafile    = 'hc,rfDC';
SUBJ(18).denoise     = 1;
SUBJ(18).cohname     = 'Std1k@6EEG/09-03-31@1034/1/';
SUBJ(18).congrROIfoi = [40 60; 36 52];

SUBJ(19).name        = 'AMN20';
SUBJ(19).rawpath     = [praw,'Project0030/'];
SUBJ(19).scanname    = 'S1k@6Ehmt/';
SUBJ(19).sessionname = '09-04-01@1036/';
SUBJ(19).runnames    = {'1/'};
SUBJ(19).emptyrun    = {'2/'};
SUBJ(19).pathname    = [pana,'Project0030/'];
SUBJ(19).datafile    = 'hc,rfDC';
SUBJ(19).denoise     = 1;
SUBJ(19).cohname     = 'Std1k@6EEG/09-04-01@1028/1/';
SUBJ(19).ecgfile     = 'AMN20ecg-aligned.mat';
SUBJ(19).ecgcomp     = [1 2 3];
SUBJ(19).congrROIfoi = [40 60; 64 80];

SUBJ(20).name        = 'LCE09';
SUBJ(20).rawpath     = [praw,'Project0030/'];
SUBJ(20).scanname    = 'S1k@6Ehmt/';
SUBJ(20).sessionname = '09-04-01@1311/';
SUBJ(20).runnames    = {'1/'};
SUBJ(20).emptyrun    = {'2/'};
SUBJ(20).pathname    = [pana,'Project0030/'];
SUBJ(20).datafile    = 'hc,rfDC';
SUBJ(20).denoise     = 1;
SUBJ(20).cohname     = 'Std1k@6EEG/09-04-01@1304/1/';
SUBJ(20).congrROIfoi = [36 56; 68 80];

clear praw pana

for k = [1:3 5:20]
  subject = SUBJ(k);
  cd(subject.pathname);
  cd('behaviour');
  if ~exist([subject.name,'behaviour.mat'], 'file')
    [rt, trl, correct, trigsmp, rs, runnr, trlid, startfix] = analyzeRT2(subject);
    cd(subject.pathname);
    cd('behaviour');
    save([subject.name,'behaviour'], 'rt', 'trl', 'correct', 'trigsmp', 'rs', 'trlid', 'runnr', 'startfix');
  else
    load([subject.name,'behaviour']);
  end
  SUBJ(k).rt      = rt;
  SUBJ(k).trl     = trl;
  SUBJ(k).correct = correct;
  SUBJ(k).trigsmp = trigsmp;
  SUBJ(k).runnr   = runnr;
  SUBJ(k).trlid   = trlid;
  SUBJ(k).startfix = startfix;
  SUBJ(k).rs      = rs;
   
  clear rt trl correct rs trigsmp;
end

cd(pwdir);

subjlist = [1:3 5:6 8:16 18:20];
