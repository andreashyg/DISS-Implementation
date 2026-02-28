function [X] = rateqsolve(T, gammas, deltas)
% Algorithm 1 in the paper. However, does not accept any matrix as input,
% but only block triangular T. For general matrices, compute Schur or real
% Schur decomp in advance and multiply the result accordingly.
% The implementation is inspired by invrat_pow.m from https://de.mathworks.com/matlabcentral/fileexchange/74317-substitution-algorithms-for-rational-matrix-equations/files/invrat_pow.m?s_tid=prof_contriblnk
arguments (Input)
    T  % block triangular (usually 1x1 or 2x2 blocks) matrix
    gammas  % coefficients of the polynomial g. Coeff at highest order first
    deltas  % coefficients of the polynomial h. Coeff at highest order first
    % whichschur % can be either 'real' or 'complex'
end

arguments (Output)
    X  % square matrix such that p(X)q(X)^-1 is approx. T
end

%% Validate input
if ~isfloat(T) || size(T, 1) ~= size(T, 2)
    error('A must be a square floating point matrix.');
end
N = size(T, 1);
if ~isnumeric(gammas) || ~isnumeric(deltas) || ~isvector(gammas) || ~isvector(deltas)
    error('C and D must be numeric vectors.');
end
m = length(gammas) - 1;
n = length(deltas) - 1;
max_deg = max(m,n);

comp_class = class(T(1,1) * gammas(1) * deltas(1));

% Make coefficient vectors row vectors.
if (size(gammas, 1) > 1)
    gammas = gammas.';
end
if (size(deltas, 1) > 1)
    deltas = deltas.';
end
gammas = [zeros(1, max_deg - m), gammas];  % add zero coeffs if gamma is lower-degree poly than delta
deltas = [zeros(1, max_deg - n), deltas];  % add zero coeffs similarily


% modify, as we want to evaluate at g(z^2), meaning, orders 0, 1, 2, ...,
% are mapped to orders 0, 2, 4
norm_gamma = zeros(1, 2*max_deg + 2);
norm_gamma(2:2:end) = gammas; % leave a zero at the start, since
% delta describe a polynomial with degree 1 higher than gamma

norm_delta = zeros(1, 2*max_deg + 2);
norm_delta(1:2:end-1) = deltas;

% We have to flip the order of the polynomial coefficients!
% Because matlab usually uses vectors that start with the highest degree.
% But in the paper (and the following code), gammas_k is the coeff for
% degree k, so we start with the lowest degree.
% Note that we do this flip AFTER the norm_gamma was created, because that
% vector does have to be in the standard matlab big-first order, because it
% is used with the function roots directly.
gammas = fliplr(gammas);
deltas = fliplr(deltas);

%% Compute Schur decomposition.

T = cast(T, comp_class);

%% Use diagonalization for normal matrices.
if isdiag(T)
    warning("Given triang. matrix was diag, solving diag entries directly")
    fT = zeros(N,1);
    for i = 1:N
       fT(i) = log(T(i,i));
    end
    
    X = diag(fT);
 
    return
end

[nu, fb, fe, fs] = block_structure(T);

%% Compute diagonal elements.
% Compute Y_{i,i}, Y^[k]_{i,i}, P^(u)_{i,i}, and Q^(v)_{i,i} for
% u = 1 : m and v = 1 : n.
l = max(m,n);
Y = zeros(N, N, l+1, comp_class);
H = zeros(N, N, comp_class);
Q = zeros(N, N, comp_class);


for i = 1:nu
  i_cpos = fb(i):fe(i);
  if fb(i) == fe(i)  % if 1x1 block
    
    j = fb(i);
    
    Y(j,j,1) = log(T(j,j)); 
   
  % if 2x2 block, solve according to Prop. 15 in DOI: 10.1016/j.laa.2018.09.010
  else
    
    alpha = (T(fb(i),fb(i)) + T(fe(i),fe(i))) / 2;
    beta = sqrt(-(T(fb(i),fb(i)) - T(fe(i),fe(i)))^2 - 4*T(fb(i),fe(i))*T(fe(i),fb(i))) / 2;
    eigval = alpha + 1i * beta;

    f_eigval = log(eigval); 

    gam = real(f_eigval);
    del = imag(f_eigval);
  
    Y(i_cpos,i_cpos,1) = (del/beta) * T(i_cpos,i_cpos) + (gam - (alpha * del) / beta) * eye(2);
  end % of evaluation of 2x2 blocks on diagonal
  
  % calculate square of Y_ii block
  Y(i_cpos, i_cpos, 2) = Y(i_cpos, i_cpos, 1) * Y(i_cpos, i_cpos, 1);
  % note: the powers are the following: Y(:,:,1) is Y, Y(:,:,2) is Y^2, but
  % Y(:,:,k) with k>2 is Y^{2(k-1)} (e.g., Y(:,:,3) is Y^4)
  
  for k = 3:l+1  % calculate (even) powers of Y_ii
    Y(i_cpos,i_cpos,k) = Y(i_cpos,i_cpos,2) * Y(i_cpos,i_cpos,k-1);
  end
  
  H(i_cpos, i_cpos) = H(i_cpos, i_cpos) + deltas(1) * eye(fs(i));
  Q(i_cpos, i_cpos) = Q(i_cpos, i_cpos) + gammas(1) * eye(fs(i));
  for k = 1:l
    H(i_cpos, i_cpos) = H(i_cpos, i_cpos) + deltas(k+1) * Y(i_cpos, i_cpos, k+1);
    Q(i_cpos, i_cpos) = Q(i_cpos, i_cpos) + gammas(k+1) * Y(i_cpos, i_cpos, k+1);
  end
  Q(i_cpos, i_cpos) = Q(i_cpos, i_cpos) - Y(i_cpos, i_cpos, 1) * H(i_cpos, i_cpos);
end


%% Compute off-diagonal elements.

for u = 1 : nu-1
  for i = 1 : nu-u
    j = i + u;
    i_cpos = fb(i):fe(i);
    j_cpos = fb(j):fe(j);
    %% Compute Y(i,j).
    
    % Compute F(i,j)^{[k]}.
    F_ij = zeros(fs(i), fs(j), max(l-1, 1), comp_class);
    for k = 2 : l
      for tt = i+1 : j-1
        tt_cpos = fb(tt):fe(tt);
        F_ij(:,:,k-1) = F_ij(:,:, k-1) +...
            Y(i_cpos, tt_cpos, 2) * Y(tt_cpos, j_cpos, k);
      end
    end
    
    phitild_ij = zeros(fs(i), fs(j), comp_class);
    for tt = i+1:j-1
      tt_cpos = fb(tt):fe(tt);
      phitild_ij = phitild_ij + Y(i_cpos, tt_cpos, 1) * Y(tt_cpos, j_cpos, 1);
    end
    
    phi_ij = zeros(fs(i), fs(j), max(l-1, 1), comp_class);
    phi_ij(:,:, 2-1) = F_ij(:,:, 2-1);
    
    for k = 3:l
      phi_ij(:,:, k-1) = F_ij(:,:, k-1) + Y(i_cpos, i_cpos, 2) * phi_ij(:,:, k-2);
    end
    
    Ehat_ij = zeros(fs(i) * fs(j), fs(i) * fs(j), l, comp_class);

    Ehat_ij(:,:, 1) = eye(fs(i) * fs(j), comp_class);

    % even powers of Ehat_ij
    for k = 2:l
      Ehat_ij(:,:, k) = kron(transpose(Y(j_cpos, j_cpos, 2)), eye(fs(i))) ...
        * Ehat_ij(:,:, k-1) ...
        + kron(eye(fs(j)), Y(i_cpos, i_cpos, k));
    end
    
    mu = zeros(fs(i), fs(i), l, comp_class);
    % calculate helper mu_k
    for k = 1:l
      mu(:,:,k) = gammas(k+1) * eye(fs(i)) + deltas(k+1) * Y(i_cpos, i_cpos, 1) ...
        + T(i_cpos, i_cpos) * (+1*deltas(k+1) * Y(i_cpos, i_cpos, 1) -1* gammas(k+1) ...
        * eye(fs(i)));
    end
  
    % calculate M_ij
    M_ij = zeros(fs(i)*fs(j), fs(i)*fs(j), comp_class);

    % first, the sum
    for k = 1:l
      M_ij(:,:) = M_ij(:,:) + kron(eye(fs(j)), mu(:,:, k)) * Ehat_ij(:,:, k);
    end
    % then, multiply
    M_ij(:,:) = M_ij(:,:) * (kron(transpose(Y(j_cpos, j_cpos, 1)),eye(fs(i))) ...
      + kron(eye(fs(j)), Y(i_cpos, i_cpos, 1)));
    % then, add last term
    M_ij(:,:) = M_ij(:,:) + kron(transpose(H(j_cpos, j_cpos)), ...
      (eye(fs(i)) + T(i_cpos, i_cpos)));

    K = zeros(fs(i), fs(j), comp_class);
    for tt = i+1:j-1
      tt_cpos = fb(tt):fe(tt);
      K = K + Y(i_cpos, tt_cpos, 1) * H(tt_cpos, j_cpos);
    end
  
    % calculate varphi_ij
    varphi_ij_helpermat = zeros(fs(i), fs(j), comp_class);
    for tt = i+1:j
      tt_cpos = fb(tt):fe(tt);
      varphi_ij_helpermat = varphi_ij_helpermat + T(i_cpos, tt_cpos) * Q(tt_cpos, j_cpos);
    end
    varphi_ij_helpermat = varphi_ij_helpermat - (eye(fs(i)) + T(i_cpos, i_cpos)) * K;

    varphi_ij = vec(varphi_ij_helpermat);
    for k=2:l
      varphi_ij = varphi_ij - kron(eye(fs(j)), mu(:,:,k)) * (Ehat_ij(:,:,k)*vec(phitild_ij) + vec(phi_ij(:,:,k-1)));
    end
    varphi_ij = varphi_ij - kron(eye(fs(j)), mu(:,:,1)) * vec(phitild_ij);
    
    % finally calculate block Y_ij
    Y(i_cpos, j_cpos, 1) = vecinv(M_ij \ varphi_ij, fs(i));

    % calculate even powers of Y_ij
    Y(i_cpos, j_cpos, 2) = Y(i_cpos, i_cpos, 1) * Y(i_cpos, j_cpos, 1) ...
      + Y(i_cpos, j_cpos, 1) * Y(j_cpos, j_cpos, 1) + phitild_ij;
    
    for k = 2:l
      Y(i_cpos, j_cpos, k+1) = Y(i_cpos, i_cpos, 2) * Y(i_cpos, j_cpos, k) ...
        + Y(i_cpos, j_cpos, 2) * Y(j_cpos, j_cpos, k) + F_ij(:,:, k-1);
    end

    % calculate H_ij
    for k = 1:l
      H(i_cpos, j_cpos) = H(i_cpos, j_cpos) + deltas(k+1) * Y(i_cpos, j_cpos, k+1);
    end

    % calculate Q_ij
    for k = 1:l
      Q(i_cpos, j_cpos) = Q(i_cpos, j_cpos) + gammas(k+1) * Y(i_cpos, j_cpos, k+1);
    end
    Q(i_cpos, j_cpos) = Q(i_cpos, j_cpos) - K ...
        - Y(i_cpos, i_cpos, 1) * H(i_cpos, j_cpos) - Y(i_cpos, j_cpos, 1) * H(j_cpos, j_cpos);

  end
end

X = Y(:,:,1);

end
