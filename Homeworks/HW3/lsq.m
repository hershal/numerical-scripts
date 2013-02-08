#!/usr/bin/octave
# Created by Hershal Bhave on 2/08/13
# For M368K HW3, § 8.1 Number 4
# Written in GNU Octave
#
# Description:
# Constructs the least squares polynomial of degree n with given
# points in the vectors x and y

function a = lsq(x, y, n)

  m = size(x);
  A = zeros(n,n);
  b = zeros(n,1);
	 
  for i=0:n
    for j=0:n
	
	A(i+1,j+1) = sum(x.^(j+i));

    endfor
    b(i+1) = sum(x(:).^(i).*y(:));

  endfor

a=A\b;
  
endfunction
