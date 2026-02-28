function [vect] = vec(mat)
%VEC transforms matrix into a column vector by stacking columns
arguments (Input)
  mat
end

arguments (Output)
  vect
end

vect = transpose(reshape(mat,1,[]));
end