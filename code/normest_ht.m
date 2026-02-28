function [est_norm] = normest_ht(X,k)
%NORMEST_HT Estimates 1-norm of X^k
arguments (Input)
    X  % square matrix
    k  % power of the matrix
end

arguments (Output)
    est_norm  % estimate of the 1-norm of X^k
end

est_norm = normest1(@Apowkmultx, 2, [], X, k);
end
