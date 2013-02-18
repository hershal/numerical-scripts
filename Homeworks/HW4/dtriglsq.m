#!/usr/bin/octave
# Created by Hershal Bhave on 2/16/13
# For M368K HW3, § 8.5 Number 7
# Written in GNU Octave
#
# Description:
# Constructs the discrete least squares trigonometric polynomial
# coefficients of degree n with the given function f (in terms of x)
# and lower and upper bounds a and b. These are nominallly -pi to pi,
# so don't count on this algorithm working for any other bounds. After
# you get the coefficients, its up to you to figure out what to do
# with them (Fourier Series is a good buzz-word suggestion). Also, if
# you provide a bad "m" value, it will automatically determine an
# appropriate m value to use in the compuatation.

function [A,B] = dtriglsq(f, n, m, a, b)

  # We're using 2*m data points to represent our polynomial.
  # Mathematically, these indices are zero-indexed since there are
  # 2*m data points in the sequence [0:2*m-1]. Remember that Octave
  # uses one-indexed vectors, so whenever I write these algorithms,
  # I use the mathematical style until I have to put them in Octave
  # vectors, in which case I usually write something like x(i+1)
  if(m<1)
    m = 2*n-1;
  endif

  A = zeros(n+1,1);
  B = zeros(n+1,1);

  x = (a:(b-a)/(2*m):b)';

  # We don't consider the last point since we want an even number of
  # data points for our algorithm. 
  x(length(x)) = [];

  y = (f(x));

  # We're using this loop for both A and B, but we only consider the
  # values in B from 1:n-1.
  for k=0:n
    A(k+1) = 1/m*sum(y(:).*cos(k*x(:)));
    B(k+1) = 1/m*sum(y(:).*sin(k*x(:)));
  endfor

endfunction
