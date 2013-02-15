#!/usr/bin/octave
# Created by Hershal Bhave on 2/08/13
# For M368K HW3, § 8.2 Number 4
# Written in GNU Octave
#
# Description:
# Constructs the continuous least squares polynomial coefficients of
# degree n with the given function f (in terms of x) and lower and
# upper bounds a and b.

function c = contlsq(f, n, a, b)

  tol=10^-6;

  i=0; j=0;

  A = zeros(n,n);
  B = zeros(n,1);

  for i=0:n
    for j=0:n

      ## !! WARNING: !!
      ## This is nasty. Consider moving this segment outside of the
      ## loop if possible. Its not good to have difficult instructions
      ## inside of a deep for-loop, especially in interpreted scripts.
      g=@(x) (x).^(i+j);
      h=@(x) (x.^i).*f(x);

      A(i+1,j+1) = quadcc(g,a,b,tol);

    endfor
    B(i+1) = quadcc(h,a,b,tol);

  endfor

c=A\B'
  
endfunction
