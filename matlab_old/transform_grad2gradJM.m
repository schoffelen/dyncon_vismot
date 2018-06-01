function [M] = transform_grad2grad(gradfrom, gradto)

if senstype(gradfrom) ~= senstype(gradto), 
  error('only transformation between same systems is possible');
end

% to construct the average location of the MEG sensors, 4 channels are needed that should  be sufficiently far apart
switch senstype(gradfrom)
case {'ctf151' 'ctf275'}
  labC = 'MZC01';
  labF = 'MZF03';
  labL = 'MLC21';
  labR = 'MRC21';
case {'ctf151_planar' 'ctf275_planar'}
  labC = 'MZC01_dH';
  labF = 'MZF03_dH';
  labL = 'MLC21_dH';
  labR = 'MRC21_dH';
case {'bti148'}
  labC = 'A14';
  labF = 'A2';
  labL = 'A15';
  labR = 'A29';
case {'bti248'}
  labC = 'A19';
  labF = 'A2';
  labL = 'A44';
  labR = 'A54';
otherwise
  % this could in principle be added to the cfg, but better is to have a more exhaustive list here
  error('unsupported MEG system for realigning, please ask on the mailing list');
end

% determine the 4 ref sensors for this individual template helmet 
indxC = strmatch(labC, gradfrom.label, 'exact'); 
indxF = strmatch(labF, gradfrom.label, 'exact'); 
indxL = strmatch(labL, gradfrom.label, 'exact'); 
indxR = strmatch(labR, gradfrom.label, 'exact'); 
if isempty(indxC) || isempty(indxF) || isempty(indxL) || isempty(indxR)
  error('not all 4 sensors were found that are needed to rotate/translate');
end
meanC = gradfrom.pnt(indxC,:);
meanF = gradfrom.pnt(indxF,:);
meanL = gradfrom.pnt(indxL,:);
meanR = gradfrom.pnt(indxR,:);

% construct two direction vectors that define the helmet orientation
dirCF = (meanF - meanC);
dirRL = (meanL - meanR);
% construct three orthonormal direction vectors
dirX = normalize(dirCF);
dirY = normalize(dirRL - dot(dirRL, dirX) * dirX);
dirZ = cross(dirX, dirY);
tra  = fixedbody(meanC, dirX, dirY, dirZ);

% determine the 4 ref sensors for the helmet that belongs to this dataset
indxC = strmatch(labC, gradto.label, 'exact'); 
indxF = strmatch(labF, gradto.label, 'exact'); 
indxL = strmatch(labL, gradto.label, 'exact'); 
indxR = strmatch(labR, gradto.label, 'exact'); 
if isempty(indxC) || isempty(indxF) || isempty(indxL) || isempty(indxR)
  error('not all 4 sensors were found that are needed to rotate/translate');
end

% construct two direction vectors that define the helmet orientation
origdirCF = gradto.pnt(indxF,:) - gradto.pnt(indxC,:);
origdirRL = gradto.pnt(indxL,:) - gradto.pnt(indxR,:);
% construct three orthonormal direction vectors
origdirX = normalize(origdirCF);
origdirY = normalize(origdirRL - dot(origdirRL, origdirX) * origdirX);
origdirZ = cross(origdirX, origdirY);
origtra = fixedbody(gradto.pnt(indxC,:), origdirX, origdirY, origdirZ);

M = inv(origtra)*tra;

function [h] = fixedbody(center, dirx, diry, dirz);
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
tra = eye(4);
tra(1:4,4)   = [-center 1]';
h = rot * tra;

function [v] = normalize(v);
v = v / sqrt(v * v');
