#!/usr/bin/octave
# Created by Hershal Bhave on 04/04/13
# For M368K HW8, § 11.1 Number 2b
# Written in GNU Octave
#
# Description: Uses the Shooting method to approximate the solution to
# a Boundary Value Problem in the form of -y''+p(x)y'+q(x)y+r(x)=0,
# given the p(x), q(x), endpoints a, b, boundary conditions aa, bb and
# the number of subintervals n
#

function [x,u,v,w,W] = shooting(p, q, r, a, b, aa, bb, n)

  h=(b-a)/n;

  u = zeros(2, n);
  v = zeros(2, n);
  k = zeros(4, n);
  kt = zeros(4, n);

  ## Step 1
  u(1, 1) = aa;
  u(2, 1) = 0;
  v(1, 1) = 0;
  v(2, 1) = 1;

  ## Step 2
  for i=1:n

    ## Step 3
    x = a+i*h;

    ## Step 4
    k(1,1) = h*u(2,i);
    k(1,2) = h*(p(x)*u(2,i)+q(x)*u(1,i)+r(x));

    k(2,1) = h*(u(2,i)+k(1,2)/2);
    k(2,2) = h*(p(x+h/2)*(u(2,i)+k(1,2)/2)+q(x+h/2)*(u(1,i)+k(1,1)/2)+r(x+h/2));

    k(3,1) = h*(u(2,i)+k(2,2)/2);
    k(3,2) = h*(p(x+h/2)*(u(2,i)+k(2,2)/2)+q(x+h/2)*(u(1,i)+k(2,1)/2)+r(x+h/2));

    k(4,1) = h*(u(2,i)+k(3,2));
    k(4,2) = h*(p(x+h)*(u(2,i)+k(3,2))+q(x+h)*(u(1,i)+k(3,1))+r(x+h));

    u(1,i+1) = u(1,i) + (k(1,1)+2*k(2,1)+2*k(3,1)+k(4,1))/6;
    u(2,i+1) = u(2,i) + (k(1,2)+2*k(2,2)+2*k(3,2)+k(4,2))/6;


    kt(1,1) = h*v(2,i);
    kt(1,2) = h*(p(x)*v(2,i)+q(x)*v(1,i));
    
    kt(2,1) = h*(v(2,i)+kt(1,2)/2);
    kt(2,2) = h*(p(x+h/2)*(v(2,i)+kt(1,2)/2)+q(x+h/2)*(v(1,i)+kt(1,1)/2));
    
    kt(3,1) = h*(v(2,i)+kt(2,2)/2);
    kt(3,2) = h*(p(x+h/2)*(v(2,i)+kt(2,2)/2)+q(x+h/2)*(v(1,i)+kt(2,1)/2));
    
    kt(4,1) = h*(v(2,i)+kt(3,2));
    kt(4,2) = h*(p(x+h)*(v(2,i)+kt(3,2))+q(x+h)*(v(1,i)+kt(3,1)));
    
    v(1,i+1) = v(1,i) + (kt(1,1)+2*kt(2,1)+2*kt(3,1)+kt(4,1))/6;
    v(2,i+1) = v(2,i) + (kt(1,2)+2*kt(2,2)+2*kt(3,2)+kt(4,2))/6;

  endfor
  
  ## Step 5
  w(1,1) = aa;
  w(2,1) = (bb-u(1,n))/(v(1,n));

  ## Step 6
  W(1,:) = u(1,1:n)+w(2,1).*v(1,1:n);
  W(2,:) = u(2,1:n)+w(2,1).*v(2,1:n);
  x = a+(0:n).*h;

endfunction
