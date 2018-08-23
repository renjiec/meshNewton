function [A] = adjacent(t)
  i = reshape(t, numel(t),1);
  j = reshape(t(:,[2 3 1]),numel(t),1);
  A = sparse(i,j,1);
end