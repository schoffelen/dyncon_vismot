function [source, parcellation] = vismot_bf_lcmv_parcellate(sourcein, tlck, varargin)

    % use an existing parcellation, but still do an svd on the projected
    % power
    parcellation = ft_getopt(varargin, 'parcellation');
    parcelparam  = ft_getopt(varargin, 'parcellationparam', 'parcellation');
    
    Nparcel = numel(parcellation.([parcelparam,'label']));
    filter  = cell(Nparcel,1);
    for k = 1:Nparcel
      sel = find(parcellation.(parcelparam)==k);
      F   = cat(1,sourcein.avg.filter{sel});
      [u,s,v] = svd(F*tlck.cov*F');
      filter{k} = u'*F;
      S{k}      = diag(s);
      U{k}      = u;
    end
    
    source = rmfield(sourcein, 'avg');
    source.parcellation = parcellation.(parcelparam);
    source.parcellationlabel = parcellation.([parcelparam,'label']);
    
    clear parcellation;
    
    parcellation.label  = source.parcellationlabel;
    parcellation.filter = filter;
    parcellation.s      = S;
    parcellation.u      = U;
    