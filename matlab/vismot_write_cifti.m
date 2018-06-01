function vismot_write_cifti(filename, data, varargin)

parameter    = ft_getopt(varargin, 'parameter',    'avg');
parcellation = ft_getopt(varargin, 'parcellation');
parcelparam  = ft_getopt(varargin, 'parcelparam', 'parcellation');

tmp = load('atlas_subparc_8k');
sourcemodel = tmp.atlas;
n = size(sourcemodel.pos,1);

% create a parcellation-structure that ft_write_cifti understands
if ~isfield(data, 'brainordinate') && ~isempty(parcellation),
  fprintf('creating a ''brainordinate'' substructure from the parcellation\n');
  if isfield(sourcemodel, 'pnt')
    data.brainordinate.pos = sourcemodel.pnt;
  elseif isfield(sourcemodel, 'pos')
    data.brainordinate.pos = sourcemodel.pos;
  end
  data.brainordinate.tri = sourcemodel.tri;
  data.brainordinate.brainstructure      = [ones(1,n/2) ones(1,n/2)*2]; % hard coded!
  data.brainordinate.brainstructurelabel = {'CORTEX_LEFT' 'CORTEX_RIGHT'};
  data.brainordinate.parcellation        = parcellation.(parcelparam)';
  data.brainordinate.parcellationlabel   = parcellation.([parcelparam,'label'])';
  
  cfg           = [];
  cfg.filetype  = 'cifti';
  cfg.parameter = parameter;
  cfg.filename  = filename;
  ft_sourcewrite(cfg, data);

elseif isfield(data, 'brainordinate')
  % ok
  ft_write_cifti(filename, data, 'parameter', parameter);
  
else
  error('you should either supply data with a ''brainordinate'' substructure, or a definition of a parcellation');
end
  