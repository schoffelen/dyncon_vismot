function dat = lrflip(dat)
% this function flips channel level data over the midline

load('/project/3011085.03/scripts/dyncon_vismot/matlab_old/lrplist.mat', 'channelcmb');

for k=1:numel(dat.label)
    idx = match_str(channelcmb(:,1), dat.label(k));
    dat.label(k) = channelcmb(idx,2);
end

