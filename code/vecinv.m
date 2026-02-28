function [mat] = vecinv(vect, colsize)
%VECINV transforms a column vector into a matrix, by slicing the vector
%into parts of fixed length 'colsize' and using those as columns
arguments (Input)
  vect
  colsize
end

arguments (Output)
  mat
end

mat = reshape(vect, [colsize, length(vect) / colsize]);
end