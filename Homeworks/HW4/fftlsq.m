#!/usr/bin/octave
# Created by Hershal Bhave on 2/16/13
# For M368K HW3, § 8.6 Number 1,2
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

function [A,B] = fftlsq(f, n, m, a, b)

  # Number of levels to use depends on the number of points we have to
  # represent.
  L = floor(log2(m));

  # Populate the points in between the two bounds.
  x = (a:(b-a)/(2*m):b)';

  # We don't consider the last point since we want an even number of
  # data points for our algorithm.
  x(length(x)) = [];

  y = (f(x));

  # Initial setup
  c=y;
  mu=exp((pi*sqrt(-1))/m);

  # From now on, C will be the new partitioned matrix and c will be
  # the old one.
  for k=0:L

      C=[
      D=eye(2^L).*mu.^(0:2^L-1);



  endfor
endfunction
