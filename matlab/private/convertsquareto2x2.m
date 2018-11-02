function [output] = convertsquareto2x2(input, colindx, input_diag)

% CONVERTSQUARETO2X2 takes an NxN matrix, and converts it into a 2x2x(M) matrix
% so that each slice is a pairwise 2x2 combination, and M represent the
% combinations as the ordered lower triangular part of the NxN matrix.
% If a second (scalar) input is provided, a 2x2x(N) matrix is returned,
% where the 2x2 combinations are returned, for one of the columns of the
% input matrix.

[m,n] = size(input);
if m~=n, error('input matrix should be square'); end
if nargin<3
  input_diag = diag(input);
end

if nargin==1
  
  tmp   = 1:n;
  indx1 = tmp(ones(n,1),:);
  indx4 = indx1';
  
  indx1 = tril(indx1,-1); % row    indices
  indx4 = tril(indx4,-1); % column indices
  
  indx1 = indx1(indx1~=0);
  indx4 = indx4(indx4~=0);
  
  % create the vector of indices
  tmp   = reshape(1:(n^2), [n n]);
  D     = diag(tmp); clear tmp;
  
  % compute only once (can become lengthy with many pairs)
  D1    = D(indx1);
  D4    = D(indx4);
  
  % allocate output
  % output = zeros(numel(indx1),2,2);
  % output(:,1,1) = input(D1); % corresponding to diagonal entries
  % output(:,2,2) = input(D4); % corresponding to diagonal entries
  % output(:,1,2) = input(D4+indx1-indx4);
  % output(:,2,1) = input(D1+indx4-indx1);
  %
  % output = permute(output, [2 3 1]);
  
  output = zeros(2,2,numel(indx1));
  output(1,1,:) = input(D1);
  output(2,2,:) = input(D4);
  output(1,2,:) = input(D4+indx1-indx4);
  output(2,1,:) = input(D1+indx4-indx1);
  
  % output = zeros(4,numel(indx1));
  % for k = 1:numel(indx1)
  %   i1 = indx1(k);
  %   i2 = indx4(k);
  %   Di1 = D(i1);
  %   Di2 = D(i2);
  %   output(:,k) = input([Di1 Di2 Di2+i1-i2 Di1+i2-i1]);
  % end
  % output = reshape(output,[2 2 numel(indx1)]);

else
  
  if numel(colindx)>1
    error('only a single column index can be provided');
  end
  output = zeros(2,2,m);
  output(1,1,:) = input_diag;%diag(input);
  output(2,2,:) = input(colindx,colindx);
  output(2,1,:) = input(colindx,:);
  output(1,2,:) = input(:,colindx);
  
end