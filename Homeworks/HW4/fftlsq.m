#!/usr/bin/octave
# Created by Hershal Bhave on 2/16/13
# For M368K HW3, § 8.6 Number 2b
# Written in GNU Octave
#
# Description:
# Constructs the least squares trigonometric polynomial coefficients
# of degree n with the given function f (in terms of x) and lower and
# upper bounds a and b using the Fast Fourier Transform Algorithm. The
# bounds are nominallly -pi to pi, so don't count on this algorithm
# working for any other bounds. After you get the coefficients, its up
# to you to figure out what to do with them (Fourier Series is a good
# buzz-word suggestion). Also, if you provide a bad "m" value, it will
# automatically determine an appropriate m value to use in the
# compuatation.

function [a,b] = fftlsq(f, n, a, b)


  A = zeros(2*n, 2*n);
  
  # Populate the points in between the two bounds.
  x = (a:(b-a)/(2*n):b)';

  # We don't consider the last point since we want an even number of
  # data points for our algorithm.
  x(length(x)) = [];

  y = (f(x));

  for i=0:2*n-1
    for j=0:2*n-1
	A(i+1,j+1) = exp(pi*1.0i/n)^(i*j);
    endfor
  endfor

  c=A*y;

  g=((exp(1.0i.*(0:2*n-1).*pi)')./n);
  a=real(g.*c);
  b=imag(g.*c);

  b(1)=0;
  b(length(b))=0;

endfunction
