function dataout = hemiflip(datain, parameter)
% flips values over the hemisphere to treat as if it is the reverse. Can for
% example be helpful to treat right hand responses as left hand responses.
% datain should be source level fieldtrip structure. parameter the
% parameter(s) to be flipped. 
% NOTE: only works when parameter is 2D array.

dataout = datain;
for k=1:numel(parameter)
    p = parameter{k};
    
    n = size(datain.(p),1)./2;
    
    dataout.(p) = datain.(p)([n+(1:n) 1:n],1);
end

            
