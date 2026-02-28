function [cs] = get_cs(m)
%GET_CS Get coefficients c_k of the pade approc of exp
%   See p. 2832
arguments (Input)
  m
end

arguments (Output)
  cs
end
% table = [ 0.5 1.]
cs = zeros(1, m+1);
for k = 0:m
    cs(k+1) = nchoosek(m, k)*(factorial(m + m - k))/(factorial(m + m));
end
cs = fliplr(cs);  % coeff for highest order first

end