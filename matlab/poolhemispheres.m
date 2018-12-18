function dataout = poolhemispheres(datain, parameter, basehemi, hemisign)
% pool the left and right hemispheres of source data. basehemi ('left', or
% 'right' determines whether the .pos info of the left or right hemisphere
% will be used. hemisign [+/-1 +/-1] determines how the data will be
% combined. NOTE: only works for structures with pos info, where -x is
% considered left hemisphere, and +x right hemisphere
% FIXME: unclear how to combine pos and dim information
error('not yet implemented')
if ~isfield(datain, 'pos')
    error('poolhemispheres only works on source structures with pos info')
end

% find rows that belong to left and right hemispheres, respectively.
left = find(datain.pos(:,1)<0);
right = find(datain.pos(:,1)>0);

dataout=datain;
% first pool parameter data
dataout.(parameter)(left,:) = hemisign(1)*dataout.(parameter)(left,:);
dataout.(parameter)(right,:) = hemisign(2)*dataout.(parameter)(right,:);

% now change pos information (test whether source structure is still
% functional if it has multiples of the same pos).
if strcmp(basehemi, 'left')
    dataout.pos(right,1) = -dataout.pos(right,1);
elseif strcmp(basehemi, 'right')
    dataout.pos(left,1) = -dataout.pos(left,1);
end

end
