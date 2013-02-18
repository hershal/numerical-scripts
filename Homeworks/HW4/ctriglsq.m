#!/usr/bin/octave
# Created by Hershal Bhave on 2/14/13
# For M368K HW3, § 8.5 Number 2
# Written in GNU Octave
#
# Description:
# Constructs the continuous least squares trigonometric polynomial
# coefficients of degree n with the given function f (in terms of x)
# and lower and upper bounds a and b. These are nominallly -pi to pi,
# so don't count on this algorithm working for any other bounds. After
# you get the coefficients, its up to you to figure out what to do
# with them (Fourier Series is a good buzz-word suggestion).

function [A,B] = ctriglsq(f, n, a, b)

  tol=10^-6;

  i=0; j=0;

  A = zeros(n,1);
  B = zeros(n,1);

  for i=0:n

    ## !! WARNING: !!
    ## This is nasty. Consider moving this segment outside of the
    ## loop if possible. Its not good to have difficult instructions
    ## inside of a deep for-loop, especially in interpreted scripts.
    g=@(x) (1/pi).*(f(x).*cos(i.*x));
    h=@(x) (1/pi).*(f(x).*sin(i.*x));
	    
    A(i+1) = quadcc(g,a,b,tol);
    B(i+1) = quadcc(h,a,b,tol);

  endfor

endfunction
