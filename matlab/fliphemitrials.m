function dataout = fliphemitrials(datain, parameter, trlinfo_column, flip_from, flip_to)

cfg=[];
cfg.trials = ismember(datain.trialinfo(:,trlinfo_column), flip_from);
fliptrials = ft_selectdata(cfg, datain);
cfg.trials = ~ismember(datain.trialinfo(:,trlinfo_column), flip_from);
nofliptrials = ft_selectdata(cfg, datain);

fliptrials = hemiflip(fliptrials, parameter);

% some hacky way to append data (doesnt work with ft_appendsource because
% of size of trialinfo)
dataout = nofliptrials;
dataout.(parameter) = [dataout.(parameter), fliptrials.(parameter)];
dataout.trialinfo = [dataout.trialinfo; fliptrials.trialinfo];

% change trialinfo
before = dataout.trialinfo(:,trlinfo_column);
after = before;
for k=1:numel(flip_from)
    after(before==flip_from(k)) = flip_to(k);
end

dataout.trialinfo(:, trlinfo_column) = after;
end