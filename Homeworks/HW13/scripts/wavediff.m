#!/usr/bin/octave
# Created by Hershal Bhave on 05/02/13
# For M368K HW13, § 12.3 Number 2
# Written in GNU Octave
#
# Description: Computes the approximate solution to the Wave Equation
# (Hyperbolic PDE), where f=u(x,0), g=(du/dt)(x,0), alpha is alpha
# in du/dt-alpha^2*(d^2u/dt^2), h is deltax, k is deltat, m is the
# length of x, and n is the length of t

function [w] = wavediff(f,g,alpha,h,k,m,n)

  i = 1:m-1;
  x = (i*h)';
  ## wold = f(x);
  ## w = wold + k*g(x);
  lambda = alpha*(k/h);

  w(:,1) = f(x);
  w(:,2) = (1-lambda^2)*f(i*h)+(lambda^2/2)*(f((i+1)*h) + \
					      f((i-1)*h)) + k*g(i*h);

  ## Not sure why this doesn't work
  A(1,1) = 2*(1-lambda^2);
  A(1,2) = lambda^2;
  
  for i=2:m-2
    A(i,i-1) = lambda^2;
    A(i,i) = 2*(1-lambda^2);
    A(i,i+1) = lambda^2;
  endfor
  
  A(m-1,m-2) = lambda^2;
  A(m-1,m-1) = 2*(1-lambda^2);
  
  for i=1:n-1
    w(:,i+2) = A*w(:,i+1) - w(:,i);;
  endfor

endfunction
