#!/usr/bin/octave
# Created by Hershal Bhave on 04/25/13
# For M368K HW12, § 12.2 Number 6
# Written in GNU Octave
#
# Description: Computes the approximate solution to the PDE using the
# Forward-Difference method, where f=u(x,0), alphasq is alpha^2 in
# du/dt-alpha^2*(d^2u/dt^2), h is deltax, k is deltat, m is the length
# of x, and n is the length of t
#

function w = forwarddiff(f,alphasq,h,k,m,n)

  i = 1:(m-1);
  x = (i*h)';
  w = f(x);

  lambda = alphasq*(k/h^2);

  A(1,1) = 1-2*lambda;
  A(1,2) = lambda;

  for i=2:m-2
    A(i,i-1) = lambda;
    A(i,i) = 1-2*lambda;
    A(i,i+1) = lambda;
  endfor

  A(m-1,m-2) = lambda;
  A(m-1,m-1) = 1-2*lambda;

  for i=1:n
    w=A*w;
  endfor

endfunction
