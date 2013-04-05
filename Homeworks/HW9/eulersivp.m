#!/usr/bin/octave
# Created by Hershal Bhave on 04/04/13
# For M368K HW8, § 11.1 Number 2a, 4a
# Written in GNU Octave
#
# Description: Uses a Second-Order Euler's Method to solve the IVP
# given in f in n steps between a and b given initial conditions on y0
# and yp0.
#

function [x,y] = eulersivp(f,n,a,b,y0,yp0)

  h = (b-a)/n;
  x = (a:h:b)';
  z(1,1) = y0;
  z(1,2) = yp0;
  
  for i=1:n
    z(i+1,1) = z(i,1) + h*z(i,2);
    z(i+1,2) = z(i,2) + h*f(x(i),z(i,1),z(i,2));
  endfor
  
  y = z(:,1);

endfunction
