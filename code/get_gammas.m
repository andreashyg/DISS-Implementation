function [gammas] = get_gammas(m)
%GET_GAMMAS Get coefficitents gamma of the function g which is used for the
% padé approximation of exp
%   In particular, gamma_k = c_2k. See p. 2839
arguments (Input)
  m
end

arguments (Output)
  gammas
end
cs = get_cs(m);  
gammas = cs(2:2:end); 
end