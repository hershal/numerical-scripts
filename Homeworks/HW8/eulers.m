#!/usr/bin/octave
# Created by Hershal Bhave on 3/25/13
# For M368K HW8, § 10.4 Number 1b
# Written in GNU Octave
#
# Description: Uses Euler's method to approximate the solution to the
# function f given an initial approximation x and number of steps n
#

function x = eulers(x, f, N)

  b = -(1/n)*f(x);
  for i=1:N
    A = jacob(x,f);
    x = x+A\b;
  endfor
endfunction
