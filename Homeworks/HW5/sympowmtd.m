#!/usr/bin/octave
# Created by Hershal Bhave on 2/23/13
# For M368K HW3, § 9.3 Number 2
# Written in GNU Octave
#
# Description: Finds the dominant eigenvalue (h) and associated
# eigenvector (x(k,:)) for the n x n matrix A and a given initial
# starting vector x using the Symmetric Power Method
#

function [x,h] = sympowmtd(A, x)

maxIter = 100;
tol = 10^-6;

k = 1;

y = x(:,k);
do

  x(:,k) = y/norm(y,2)
  y = A*x(:,k);

  h(k) = x(:,k)'*y;

  err(k) = norm(y/h(k) - x(:,k));
  
  k++;

until k>maxIter || err(k-1)<tol || abs(h(k-1))<tol

if abs(h(k-1))<tol
  error("Eigenvalue converged to zero, select a new initial guess");
endif

if k>maxIter
   error("Maximum iterations was exceeded, select a new initial guess");
endif

endfunction
