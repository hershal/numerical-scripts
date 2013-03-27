#!/usr/bin/octave
# Created by Hershal Bhave on 3/25/13
# For M368K HW8, § 10.5 Number 1b
# Written in GNU Octave
#
# Description: Uses the Runge-Kutta Method of order 4 to approximate
# the solution to the function f given an initial approximation x
#

function x = rungekutta(x, f, N)

  h = 1/N;
  b = -h*f(x);

  if size(x,1) < size(x,2)
    x=x';
  endif

  for i=1:N
      A = jacob(x,f);
      k(:,1) = A\b;

      A = jacob(x+1/2.*k(:,1),f);
      k(:,2) = A\b;

      A = jacob(x+1/2.*k(:,2),f);
      k(:,3) = A\b;

      A = jacob(x+k(:,3),f);
      k(:,4) = A\b;

      x = x + (k(:,1)+2*k(:,2)+2*k(:,3)+k(:,4))/6;
  endfor

endfunction
