#!/usr/bin/octave
# Created by Hershal Bhave on 04/25/13
# For M368K HW12, § 12.2 Number 10
# Written in GNU Octave
#
# Description: Uses the Finite-Difference algorithm to approximate the
# solution to a given Poisson equation f, its boundary conditions g, x
# boundaries a and b, y boundaries c and d, and grid dimensions m and n
#

function [w, u] = crank(l,T,alpha,m,N)

  ## Step 1
  h = l/m;
  k = T/N;
  lambda = alpha^2*k/h^2;
  w(m) = 0;

  ## Step 2
  for i=1:m-1
    w(i) =f(i*h);
  endfor

  ## Step 3
  ## Solving the linear system
  l(1) = 1 + lambda;
  u(1) = -lambda/(2*l(1));

  ## Step 4
  for i=2:m-2
      l(i) = 1 + lambda + lambda*u(i-1)/2;
      u(i) = -lambda/(2*l(i));
  endfor
  
  ## Step 5
  l(m-1) = 1 + lambda + lambda*u(m-2)/2;

  ## Step 6
  for j=1:N

    ## Step 7
    t = j*k;
    z(1) = ((1-lambda)*w(1) + (lambda/2)*w(2))/l(1);

    ## Step 8
    for i=2:m-1
      z(i) = ((1-lambda)*w(i) + (lambda/2)*(w(i+1) + w(i-1) + z(i-1)))/l(i);
    endfor

    ## Step 9
    w(m-1) = z(m-1);

    ## Step 10
    for i=m-2:-1:1
	w(i) = z(i) - u(i)*w(i+1);
    endfor
    
    ## Step 11
    for i=1:m-1
      x = i*h;
    endfor

  endfor

endfunction
