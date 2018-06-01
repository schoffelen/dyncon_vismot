function [data1,data2,data3,data4,data5] = denoise_reject(subject,varargin)

%collects data in the 5 conditions for a given subject
%estimating digital weights on selected trials across all conditions
%or using the weights specified in the other input arguments
%
%use as:
% [data1,data2,data3,data4,data5] = denoise_reject(subject)
%
%denoising should be done across all conditions to make life easier
%otherwise a condition specific balancing will be introduced

if length(varargin)>0,
  for k = 1:length(varargin)
    tmp  = inputname(1+k);
     
    data = varargin{k};
    trl  = findcfg(data.cfg, 'trl');
    w    = findcfg(data.cfg, 'pca');
  
    [data,refdata] = doPreprocessingNew(subject, trl);

    cfg     = [];
    cfg.pca = w;
    data    = denoise_pca(cfg, data, refdata);
    clear refdata

    cfg     = [];
    cfg.channel = channelselection(varargin{k}.label, data.label);
    data    = preprocessing(cfg, data);

    cfg            = [];
    cfg.resamplefs = 256;
    cfg.blc        = 'yes';
    cfg.blcwindow  = [-0.5 0];
    cfg.detrend    = 'no';
    data           = resampledata(cfg, data);

    %append headposition data if applicable
    [a,b] = match_str(data.label, varargin{k}.label);
    cfg   = [];
    cfg.channel = varargin{k}.label(setdiff(1:length(varargin{k}.label), b(:)'));
    tmpdata     = struct2double(preprocessing(cfg, varargin{k}));
    tmpdata.time = data.time;
    data        = appenddata([], data, tmpdata);
    
    if strcmp(tmp,'data1')
      data1 = struct2single(data);
    elseif strcmp(tmp,'data2')
      data2 = struct2single(data);
    elseif strcmp(tmp,'data3')
      data3 = struct2single(data);
    elseif strcmp(tmp,'data4')
      data4 = struct2single(data);
    elseif strcmp(tmp,'data5')
      data5 = struct2single(data);
    end
    clear data
  end  
  if ~exist('data1', 'var'), data1 = []; end
  if ~exist('data2', 'var'), data2 = []; end
  if ~exist('data3', 'var'), data3 = []; end
  if ~exist('data4', 'var'), data4 = []; end
  if ~exist('data5', 'var'), data5 = []; end
else
  %this is the original code, but as it turns out, line-noise has not been adequately removed
  %part of the preprocessing has to be done; i.e. application of a proper dftfilter to the trials which are to be
  %kept eventually. the weights do not have to be recomputed, nor does the rejectvisual have to be done. the stuff
  %before the else tries to achieve this, by loading in the data, and reconstructing it from the trl per condition,
  %do pca projecting according to the weights in the data, and dftfiltering it, before downsampling

  refdata = collectRefdata(subject);
  megg    = channelselection({'MEGREFG'}, refdata.label);
  mega    = channelselection({'MEGREFA'}, refdata.label);
  megr    = channelselection({'MEGREF' '-MEGREFA' '-MEGREFG'}, refdata.label);
  
  %select a subset of trials which will be used for the weights computation
  cfg         = [];
  cfg.method  = 'summary';
  cfg.metric  = 'var';
  cfg.channel = megg;
  cfg.keepchannel = 'yes';
  tmp         = rejectvisual(cfg, refdata);
  close all;
  cfg.metric  = '1/var';
  tmp         = rejectvisual(cfg, tmp);
  close all;
  trl1        = findcfg(tmp.cfg, 'trl');
  cfg.channel = mega;
  cfg.metric  = 'var';
  tmp         = rejectvisual(cfg, refdata);
  close all;
  cfg.metric  = '1/var';
  tmp         = rejectvisual(cfg, tmp);
  close all;
  trl2        = findcfg(tmp.cfg, 'trl');
  cfg.channel = megr;
  cfg.metric  = 'var';
  tmp         = rejectvisual(cfg, refdata);
  close all;
  cfg.metric  = '1/var';
  tmp         = rejectvisual(cfg, tmp);
  close all;
  trl3        = findcfg(tmp.cfg, 'trl');
  clear tmp
  
  trl           = findcfg(refdata.cfg, 'trl');
  [trl2,ia,ib]  = intersect(trl,intersect(intersect(trl1,trl2,'rows'),trl3,'rows'),'rows');
  
  %collect the MEG data and put trials in the correct order
  [data,newtrl] = collectData(subject, trl2);
  [trl3,ia,ib]  = intersect(trl,newtrl,'rows');
  refdata.trial = refdata.trial(ia);
  refdata.time  = refdata.time(ia);
  data.trial    = data.trial(ib);
  data.time     = data.time(ib);
  
  %compute the weights 
  cfg            = [];
  cfg.refchannel = channelselection({'MEGREF'}, refdata.label);
  cfg.channel    = channelselection({'MEG'}, data.label);
  cfg.zscore     = 'yes';
  cfg.truncate   = 1e-8;
  data           = denoise_pca(cfg, data, refdata);
  w              = data.cfg.pca;
  
  %now collect all trials again and do pca-subtraction
  cfg     = [];
  cfg.pca = w;
  
  cfg2         = [];
  cfg2.method  = 'summary';
  cfg2.channel = 'MEG';
  
  cfg3            = [];
  cfg3.resamplefs = 256;
  cfg3.blc        = 'yes';
  cfg3.blcwindow  = [-0.5 0];
  cfg3.detrend    = 'no';
  
  [data,refdata] = doPreprocessing(subject,1);
  data           = denoise_pca(cfg,data,refdata);
  data           = rejectvisual(cfg2, data);
  close all;
  data1          = resampledata(cfg3, data);
  data1          = struct2single(data1);
  clear data refdata;
  [data,refdata] = doPreprocessing(subject,2);
  data           = denoise_pca(cfg,data,refdata);
  data           = rejectvisual(cfg2, data);
  close all;
  data2          = resampledata(cfg3, data);
  data2          = struct2single(data2);
  clear data refdata;
  [data,refdata] = doPreprocessing(subject,3);
  data           = denoise_pca(cfg,data,refdata);
  data           = rejectvisual(cfg2, data);
  close all;
  data3          = resampledata(cfg3, data);
  data3          = struct2single(data3);
  clear data refdata;
  [data,refdata] = doPreprocessing(subject,4);
  data           = denoise_pca(cfg,data,refdata);
  data           = rejectvisual(cfg2, data);
  close all;
  data4          = resampledata(cfg3, data);
  data4          = struct2single(data4);
  clear data refdata;
  [data,refdata] = doPreprocessing(subject,5);
  data           = denoise_pca(cfg,data,refdata);
  data           = rejectvisual(cfg2, data);
  close all;
  data5          = resampledata(cfg3, data);
  data5          = struct2single(data5);
  clear data refdata;
end



%data           = appenddata([], data1, data2, data3, data4, data5);
%
%data           = struct2double(data);
%cfg            = [];
%%cfg.runica.pca = 100;
%cfg.method     = 'runica';
%comp           =  componentanalysis(cfg, data);

%to do save comp along with data
