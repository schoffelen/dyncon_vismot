function dataout = poolhemispheres(datain, parameter, basehemi, hemisign, template)
% pool the left and right hemispheres of source data. basehemi ('left', or
% 'right' determines whether the .pos info of the left or right hemisphere
% will be used. hemisign [+/-1 +/-1] determines how the data will be
% combined. NOTE: only works for structures with pos info, where -x is
% considered left hemisphere, and +x right hemisphere
% FIXME: unclear how to combine pos and dim information

if ~isfield(datain, 'pos')
    error('poolhemispheres only works on source structures with pos info')
end

template = ft_convert_units(template, 'mm');

% find rows that belong to left and right hemispheres, respectively.
left = find(template.pos(:,1)<0);
leftpos = template.pos(left,:);
right = find(template.pos(:,1)>0);
rightpos = template.pos(right,:);

% order the right hemisphere positions according to the mirror position in
% the left hemisphere
[~, tmpidx] = ismember(rightpos, [-1 1 1].*leftpos, 'rows');
right = right(tmpidx);

% take the mean over mirror locations of the left and right hemisphere
dataout=datain;
dataout.(parameter)(left,:) = (hemisign(1)*datain.(parameter)(left,:) + hemisign(2)*datain.(parameter)(right,:))/2;
dataout.(parameter)(right,:) = dataout.(parameter)(left,:);

% convert one hemisphere's inside to zeros.
dataout = datain;
if strcmp(basehemi, 'left')
    dataout.inside(right)=0;
elseif strcmp(basehmi, 'right')
    dataout.inside(left)=0;
end

end
