function index = find_dipoleindex(data, location)
% find the dipole index number from a .pos field of a fieldtrip source data
% structure, matching a given location. location has to be in mm.

% make sure units are the same
data = ft_convert_units(data, 'mm');

% make sure there are no rounding inconsistencies
pos = round(data.pos);
location = round(location); 

index = find(ismember(pos, location, 'rows'));
end



