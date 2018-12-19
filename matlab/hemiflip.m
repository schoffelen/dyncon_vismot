function dataout = hemiflip(datain, parameter)
% flips values over the hemisphere to treat as if it is the reverse. Can for
% example be helpful to treat right hand responses as left hand responses.
% datain should be source level fieldtrip structure. parameter the
% parameter(s) to be flipped.
% NOTE: only works when parameter is 2D array.

ntrials = size(datain.(parameter),2);

dataout = datain;
if iscell(parameter)
    for k=1:numel(parameter)
        p = parameter{k};
        if isfield(datain, 'tri')
            n = size(datain.(p),1)./2;
            dataout.(p) = datain.(p)([n+(1:n) 1:n],:);
        elseif isfield(datain, 'dim')
            for m=1:ntrials
                dataout.(p)(:,m) = reshape(flip(reshape(datain.(p)(:,m), datain.dim),1),[],1);
            end
        end
    end
else
    if isfield(datain, 'tri')
        n = size(datain.(parameter),1)./2;
        dataout.(parameter) = datain.(parameter)([n+(1:n) 1:n],:);
    elseif isfield(datain, 'dim')
        for m=1:ntrials
            dataout.(parameter)(:,m) = reshape(flip(reshape(datain.(parameter)(:,m), datain.dim),1),[],1);
        end
    end
end

