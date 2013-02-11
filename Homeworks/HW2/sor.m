#!/usr/bin/octave
# Created by Hershal Bhave on 1/31/13
# For M368K HW2, § 7.4 Number 2
# Written in GNU Octave

function x = sor(A, b, w, x)

  i = 0;
  tol = 10^-6;
  maxIter = 100;

  do
    n = length(A);
    k = size(x,2)+1;
    x(:,k) = zeros(n,1);

    D = eye(n).*A;
    L = -tril(A, -1);
    U = -triu(A,  1);

    T = inverse(D-w*L)*((1-w)*D + w*U);
    c = w*inverse(D-w*L)*b;

    # The magic happens here
    x(:,k) = T*x(:,k-1) + c;

    i++;
  until((abs(max(x(:,k)-x(:,k-1)))<tol) || (i>maxIter))

  if (i>maxIter)
    fprintf("SOR ERROR: Maximum iterations exceeded. Consider decreasing tolerance increasing maximum iterations, or tweak your input parameters\n");
  endif

endfunction

