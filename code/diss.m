function logA = diss(A, whichschur)
%DISS Algorithm 2 in the paper
arguments (Input)
    A  % square matrix
    whichschur  % 'complex' or 'real'
end

arguments (Output)
    logA  % approximate matrix logarithm of the given matrix A
end

%% Define constants
% from Table 2
theta_G = dictionary(1, 1.1003511163692342e-5, ...
       2, 2.4012849957497128e-3, ...
       3, 2.7099573188927441e-2, ... 
       5, 2.6059916466908718e-1, ...
       7, 6.5282885430846634e-1, ...
       9, 9.0572865457020838e-1);

maxroots = 100;  % same value as in matlabs logm

% Compute the Schur decomposition A= UTU*
[U, T] = schur(A, whichschur);

% assert spectrum not on negative real line (code copied from matlabs logm
% implementation)
ei = ordeig(T);
warns = any(ei == 0);
if any(real(ei) < 0 & imag(ei) == 0 )
    warns = true;
end
if warns
    warning("given matrix has real nonpositive eigenvalues");
end

s = 0;
m = 0;
B = T;

while abs(eigs(B - eye(size(B,1)), 1, 'largestabs')) > theta_G(9) && s < maxroots
    B = sqrtm(B);
    s = s + 1;
end

if s == maxroots
    warning("Max number of roots already reached in first loop. ..." + ...
        "Will still proceed with the algorithm, but skip any further roots.")
end

%% Parameter selection loop
while m == 0
    mu_4 = nthroot(normest_ht(B - eye(size(B,1)), 4), 4);
    a_3 = max(nthroot(normest_ht(B - eye(size(B,1)), 3), 3), mu_4);
    i = 9;
    while i >= 3
        if a_3 <= theta_G(i)
            m = i;
        end
        i = i - 2;
    end
    if m == 0
        a_4 = max(mu_4, nthroot(normest_ht(B - eye(size(B,1)), 5), 5));
        zeta = min(a_3, a_4);
        if zeta <= theta_G(9)
            if zeta <= theta_G(7)
                m = 7;
            else
                m = 9;
            end
        elseif s == maxroots
            warning("Max number of square roots was not enough, will still " + ...
                "proceed to evaluate but with no further roots taken. Setting m to max, m=9")
            m = 9;

            break
        else
            B = sqrtm(B);
            s = s + 1;
        end
    end
end

logA = 2^s * U * rateqsolve(B, get_gammas(m), get_deltas(m)) * U';

end
