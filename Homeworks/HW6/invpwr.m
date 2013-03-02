#!/usr/bin/octave
# Created by Hershal Bhave on 2/28/13
# For M368K HW3, § 9.3 Number 4
# Written in GNU Octave
#
# Description: Finds the dominant eigenvalue (h) and associated
# eigenvector (x(k,:)) for the n x n matrix A and a given initial
# starting vector x using the Inverse Power Method
#

function [x,h] = invpwr(A, x, q)

maxIter = 100;
tol = 10^-6;

k=1;

[tmp,p] = max(abs(x(:,k)));
x(:,k) = x(:,k)/tmp;

do
  y = (A - (q.*eye(size(A,1))))\x(:,k);
  # Throw away tmp because we need its true value, not its abs value;
  # we only want p(k) because it contains the index of the largest
  # component of y 
  mu = y(p);
  [tmp,p] = max(abs(y));
  x(:,k+1)=y/y(p);
  x(:,k+1)

  err(k) = norm(x(:,k)-x(:,k+1),inf);
  err(k)
  h = (1/mu)+q
  k
printf("\n");
  k++;

until k>maxIter || err(k-1)<tol

if k>maxIter
   error("Maximum iterations was exceeded, select a new initial guess");
endif

h = (1/mu)+q;

endfunction
