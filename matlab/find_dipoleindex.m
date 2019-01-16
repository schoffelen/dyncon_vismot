function index = find_dipoleindex(data, location)
% find the dipole index number from a .pos field of a fieldtrip source data
% structure, matching a given location. location has to be in cm.

% transform units to mm, while assuming they are in cm.
data = ft_convert_units(data, 'mm');
location = location*10;

% make sure there are no rounding inconsistencies
pos = round(data.pos);
location = round(location); 

[~,index] = (ismember(location,pos, 'rows'));
end



