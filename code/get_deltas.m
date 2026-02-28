function [deltas] = get_deltas(m)
%GET_DELTAS Get coefficients delta of the function h which is used for the
% padé approximation of exp
%   In particular, delta_k = c_2k+1. See p. 2839
arguments (Input)
  m
end

arguments (Output)
  deltas
end
cs = get_cs(m);
deltas = cs(1:2:end-1); 
end
