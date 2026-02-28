function [out] = Apowkmultx(flag, x, mat, k)
%GET_AKXFUN Function that returns A^k*x without calculating A^k
%   To be passed to normest1
arguments (Input)
  flag  % what to output, can be 'dim', 'real', 'notransp' or 'transp'
  x  % vector to multiply A^k with
  mat  % the matrix A that is to be powered to the k and multiplied with x
  k  % the power of A
end

arguments (Output)
  out
end

if isequal(flag, 'dim')
  out = size(mat, 1);
  return
elseif isequal(flag, 'real')
  out = isreal(mat);
  return
elseif isequal(flag, 'notransp')
  out = x;
  for i = 1:k
    out = mat*out;
  end
  return
elseif isequal(flag, 'transp')
  out = x;
  for i = 1:k
    out = mat'*out;
  end
  return
else
  error('no valid flag')
end

end