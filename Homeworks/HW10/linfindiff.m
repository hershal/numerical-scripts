#!/usr/bin/octave
# Created by Hershal Bhave on 04/11/13
# For M368K HW10, § 11.3 Number 2a
# Written in GNU Octave
#
# Description: Uses the Linear Finite-Difference method to approximate
# the solution within the interval of the Boundary Value Problem in
# the form of -y''+p(x)y'+q(x)y+r(x)=0, given p(x), q(x), r(x),
# endpoints a, b, boundary conditions alpha, beta and subintervals n.
#

function [x,w] = linfindiff(p, q, r, a, b, alpha, beta, n)

  ai = bi = ci = di = zeros(n, 1);

  h = (b - a)/(n+1);
  x = a + h;
  ai(1) = 2 + (h^2)*q(x);
  bi(1) = -1 + (h/2)*p(x);
  di(1) = -h^2*r(x) + (1 + (h/2)*p(x))*alpha;

  for i=2:n-1
    x = a + i*h;
    ai(i) = 2 + (h^2)*q(x);
    bi(i) = -1 + (h/2)*p(x);
    ci(i) = -1 - (h/2)*p(x);
    di(i) = -(h^2)*r(x);
  endfor

  x = b - h;
  ai(n) = 2 + (h^2)*q(x);
  ci(n) = -1 - (h/2)*p(x);
  di(n) = -h^2*r(x) + (1 - (h/2)*p(x))*beta;

  l(1) = ai(1);
  u(1) = bi(1)/ai(1);
  z(1) = di(1)/l(1);

  for i=2:n-1
      l(i) = ai(i) - ci(i)*u(i-1);
      u(i) = bi(i)/l(i);
      z(i) = (di(i) - ci(i)*z(i-1))/l(i);
  endfor

  l(n) = ai(n) - ci(n)*u(n-1);
  z(n) = (di(n) - ci(n)*z(n-1))/l(n);
  
  w(n+1) = beta;
  w(n) = z(n);

  for i=n-1:-1:1
      w(i) = z(i)-u(i)*w(i+1);
  endfor

  x = a:h:b;
  w = [alpha w];

endfunction
