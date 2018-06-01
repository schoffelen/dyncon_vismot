function plot4D(input, param, cax)

if nargin<3
  cax = [];
end

dim = pos2dim3d(input.pos);
input.(param) = reshape(input.(param), [prod(dim) numel(input.freq)]);
inside = input.inside;
input.inside = input.(param);
input.inside(:) = 0;
input.inside(inside,:) = 1;
input.inside = input.inside>0;
input.dim = [dim numel(input.freq)];

cfg = [];
cfg.funparameter = param;
cfg.interactive = 'yes';
cfg.funcolormap = 'jet';
if ~isempty(cax)
  cfg.funcolorlim = cax;
end
figure;ft_sourceplot(cfg, input);
