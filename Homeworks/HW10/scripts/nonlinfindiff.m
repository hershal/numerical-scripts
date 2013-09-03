#!/usr/bin/octave
# Created by Hershal Bhave on 04/11/13
# For M368K HW10, § 11.4 Number 4a
# Written in GNU Octave
#
# Description: Uses the Nonlinear Finite-Difference method to approximate
# the solution at values within the interval of the Boundary Value
# Problem in the form of -y''+p(x)y'+q(x)y+r(x)=0, given p(x), q(x),
# r(x), endpoints a, b, boundary conditions alpha, beta, the number
# of subintervals n, and maximum iterations M.
#

function [x,w] = nonlinfindiff(f, fy, fyp, a, b, alpha, beta, n, M)

  if n<2
     error("n must be >=2\n");
  endif

  tol = 10^-6;
  w = zeros(1, n+1);
  ai = bi = ci = di = zeros(1,n);
  h = (b - a)/(n+1);
  w(n+1) = beta;

  for i=1:n
    w(i) = alpha + i*((beta-alpha)/(b-a))*h;
  endfor
  
  k=1;
  do 
    x = a + h;
    t=(w(2)-alpha)/(2*h);
    ai(1) = 2 + h^2*fy(x,w(1),t);
    bi(1) = -1 + (h/2)*fyp(x,w(1),t);
    di(1) = -(2*w(1) - w(2) - alpha + h^2*f(x,w(1),t));
    
    for i=2:n-1
      x = a + i*h;
      t = (w(i+1)-w(i-1))/(2*h);
      ai(i) = 2 + h^2*fy(x,w(i),t);
      bi(i) = -1 + (h/2)*fyp(x,w(i),t);
      ci(i) = -1 - (h/2)*fyp(x,w(i),t);
      di(i) = -(2*w(i) - w(i+1) - w(i-1) + h^2*f(x,w(i),t));
    endfor

    x = b - h;
    t = (beta - w(n-1))/(2*h);
    ai(n) = 2 + h^2*fy(x,w(n),t);
    ci(n) = -1 - (h/2)*fyp(x,w(n),t);
    di(n) = -(2*w(n) - w(n-1) - beta + h^2*f(x,w(i),t));

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
    
    v(n) = z(n);
    w(n) += v(n);

    for i=n-1:-1:1
      v(i) = z(i)-u(i)*v(i+1);
      w(i) = w(i)+v(i);
    endfor

    k++;
    until k>M || norm(v) < tol

  if(k>M)
    printf("max iteration exceeded\n");
  endif

  x = a:h:b;
  w = [alpha w];

endfunction
