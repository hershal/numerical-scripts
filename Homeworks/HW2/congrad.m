#!/usr/bin/octave
# Created by Hershal Bhave on 1/31/13
# For M368K HW2, § 7.6 Number 6
# This is an incomplete implementation used only to answer the question
# Written in GNU Octave

function [x, r] = congrad(A, b, C_inv, x)

  i = 0;
  tol = 10^-6;
  maxIter = 100;

  r = b - A*x(:,1);
  w = C_inv*r;
  v = C_inv*w;

  alpha = sum(w.^2);

  do
    k = size(x,2)+1;

    u = A*v;
    t = (alpha)/(sum(v.*u));
    x(:,k) = x(:,k-1)+t*v;
    r = r - t*u;
    w = C_inv*r;
    beta = sum(w.^2);

    if ((abs(beta) < tol) && (norm(r) < tol))
      break;
    endif

    s = beta/alpha;
    v = C_inv*w + s*v;
    alpha=beta;

    i++;
  until(i>maxIter)
  
endfunction
