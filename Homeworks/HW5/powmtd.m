#!/usr/bin/octave
# Created by Hershal Bhave on 2/23/13
# For M368K HW3, § 9.3 Number 2,14
# Written in GNU Octave
#
# Description: Finds the dominant eigenvalue (h) and associated
# eigenvector (x(k,:)) for the n x n matrix A and a given initial
# starting vector x using the Power Method.
#

function [x,h] = powmtd(A, x)

maxIter = 100;
tol = 10^-6;

k=1;

# Initial scaling; we want to keep track of p
[tmp,p(k)] = max(abs(x(:,k)));
x(:,k) = x(:,k)/tmp;

do

  y = A*x(:,k);

  # Throw away tmp because we need its true value, not its abs value;
  # we only want p(k) because it contains the index of the largest
  # component of y 
  [tmp,p(k+1)] = max(abs(y));

  h(k) = y(p(k+1));
  x(:,k+1) = y/h(k);

  err(k) = norm(A*x(:,k)-h(k)*x(:,k),inf)/abs(h(k));
  
  k++;

until k>maxIter || err(k-1)<tol || abs(h(k-1))<tol

if abs(h(k-1))<tol
  error("Eigenvalue converged to zero, select a new initial guess");
endif

if k>maxIter
   error("Maximum iterations was exceeded, select a new initial guess");
endif

endfunction
