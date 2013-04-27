#!/usr/bin/octave
# Created by Hershal Bhave on 04/25/13
# For M368K HW12, § 12.2 Number 10
# Written in GNU Octave
#
# Description: Computes the approximate solution to the PDE using the
# Crank-Nicolson method, where f=u(x,0), alphasq is alpha^2 in
# du/dt-alpha^2*(d^2u/dt^2), h is deltax, k is deltat, m is the length
# of x, and n is the length of t
#

function [w,u] = crankn2(f,alphasq,h,k,m,n)

  i = 1:(m-1);
  x = (i*h)';
  w = f(x);

  lambda = alphasq*(k/h^2)

  A(1,1) = 1+lambda;
  A(1,2) = -lambda/2;
  B(1,1) = 1-lambda;
  B(1,2) = lambda/2;
  
  for i=2:m-2
    A(i,i-1) = -lambda/2;
    A(i,i) = 1+lambda;
    A(i,i+1) = -lambda/2;
    B(i,i-1) = lambda/2;
    B(i,i) = 1-lambda;
    B(i,i+1) = lambda/2;
  endfor

  A(m-1,m-2) = -lambda/2;
  A(m-1,m-1) = 1+lambda;
  B(m-1,m-2) = lambda/2;
  B(m-1,m-1) = 1-lambda;

  for i=1:n
      w=A\B*w;
  endfor

endfunction
