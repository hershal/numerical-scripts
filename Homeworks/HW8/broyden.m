#!/usr/bin/octave
# Created by Hershal Bhave on 3/25/13
# For M368K HW8, § 10.2 Number 2
# Written in GNU Octave
#
# Description: Uses the Broyden Algorithm to approximate the solution
# to the nonlinear system given an initial approximation vector x and
# an array of functions specified in a function handler.
#

function x = broyden(x, f, n)

  pkg load optim;

  tol = 10^-6;

  A = jacob(x,f);
  v = f(x);
  A = inv(A);
  s = -A*v;
  x = x+s;
  k=1;

  while(k<n && norm(s)>tol)
    w=v;
    v=f(x);
    y=v-w;
    z=-A*y;
    p=-s'*z;
    u=(s'*A)';
    A=A+(1/p)*(s+z)*u';
    s=-A*v;
    x=x+s;
    k++;
  endwhile

endfunction
