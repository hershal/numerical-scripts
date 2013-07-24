#!/usr/bin/octave
# Created by Hershal Bhave on 04/18/13
# For M368K HW11, § 12.1 Number 2
# Written in GNU Octave
#
# Description: Uses the Finite-Difference algorithm to approximate the
# solution to a given Poisson equation f, its boundary conditions g, x
# boundaries a and b, y boundaries c and d, and grid dimensions m and n
#

function [x,y,w,l] = poissonFinDiff(a,b,c,d,f,g,m,n,N)

  if (m<3 || n<3)
    error("m and n must be greater than 3");
  endif

  tol = 10^-10;

  # Step 1
  h = (b-a)/n;
  k = (d-c)/m;

  # Step 2, 3
  x = a + [1:(n-1)]*h;
  y = c + [1:(m-1)]*k;

  # Step 4
  w = zeros(n-1,m-1);

  # Step 5
  lam = h^2/k^2;
  mu = 2*(1+lam);
  l=0;

  # Step 6
  do
    # Step 7
    z = (-h^2*f(x(1),y(m-1)) + g(a,y(m-1)) + lam*g(x(1),d) + lam*w(1,m-2) + w(2,m-1))/mu;
    NORM = abs(z-w(1,m-1));
    w(1,m-1) = z;

    # Step 8
    for i = 2:(n-2)
      z = (-h^2*f(x(i),y(m-1)) + lam*g(x(i),d) + w(i-1,m-1) + w(i+1,m-1) + lam*w(i,m-2))/mu;
      NORM = max(abs(w(i,m-1)-z), NORM);
      w(i,m-1) = z;
    endfor

    # Step 9
    z = (-h^2*f(x(n-1),y(m-1)) + g(b,y(m-1)) + lam*g(x(n-1),d) + w(n-2,m-1) + lam*w(n-1,m-2))/mu;
    NORM = max(abs(w(n-1,m-1)-z), NORM);
    w(n-1,m-1) = z;

    # Step 10
    for j = (m-2):-1:2

      # Step 11
      z = (-h^2*f(x(1),y(j)) + g(a,y(j)) + lam*w(1,j+1) + lam*w(1,j-1) + w(2,j))/mu;
      NORM = max(abs(w(1,j)-z), NORM);
      w(1,j) = z;

      # Step 12
      for i=2:(n-2)
	z = (-h^2*f(x(i),y(j)) + w(i-1,j) + lam*w(i,j+1) + w(i+1,j) + lam*w(i,j-1))/mu;
	NORM = max(abs(w(i,j)-z), NORM);
	w(i,j) = z;
      endfor

      # Step 13
      z = (-h^2*f(x(n-1),y(j)) + g(b,y(j)) + w(n-2,j) + lam*w(n-1,j+1) + lam*w(n-1,j-1))/mu;
      NORM = max(abs(w(n-1,j)-z), NORM);
      w(n-1,j) = z;

    endfor

    # Step 14
    z = (-h^2*f(x(1),y(1)) + g(a,y(1)) + lam*g(x(1),c) + lam*w(1,2) + w(2,1))/mu;
    NORM = max(abs(w(1,1)-z), NORM);
    w(1,1) = z;

    # Step 15
    for i=2:(n-2)
      z = (-h^2*f(x(i),y(1)) + lam*g(x(i),c) + w(i-1,1) + lam*w(i,2) + w(i+1,1))/mu;
      NORM = max(abs(w(i,1)-z), NORM);
      w(i,1) = z;
    endfor

    # Step 16
    z = (-h^2*f(x(n-1),y(1)) + g(b,y(1)) + lam*g(x(n-1),c) + w(n-2,1) + lam*w(n-1,2))/mu;
    NORM = max(abs(w(n-1,1)-z), NORM);
    w(n-1,1) = z;

    l++;
    until l>N || NORM<tol

  if(l>N)
    error("Maximum iterations exceeded\n");
  endif
  
endfunction
