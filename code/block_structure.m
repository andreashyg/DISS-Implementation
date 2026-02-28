function [m,fb,fe,s] = block_structure(T)
% Copied from https://ch.mathworks.com/matlabcentral/fileexchange/74317-substitution-algorithms-for-rational-matrix-equations?s_tid=prof_contriblnk

% [m,fb,fe,s]=block_structure(T) computes the block structure of
% of the quasi-upper triangular matrix T
% m is the number of diagonal blocks
% fb is an array containing the begin of each block
% fe is an array containing the end of each block
% s is an array containing the sizes of the diagonal blocks
  n=length(T);
  tol=1e-15;
  i=1;
  j=1;
  while(i<n)
    v(j)=i;
    if(abs(T(i+1,i))<tol)
      i=i+1;s(j)=1;
    else
      i=i+2;s(j)=2;
    end
    j=j+1;
  end
  if(i==n)
    v(j)=n;s(j)=1;
  end
  m=length(v);
  for i=1:m
    fb(i)=v(i);
    fe(i)=v(i)+s(i)-1;
  end
end
