#!/usr/bin/octave
# Created by Hershal Bhave on 2/28/13
# For M368K HW3, § 9.3 Number 14
# Written in GNU Octave
#
# Description: Finds the deflated matrix B for the matrix using the
# largest eigenvalue's eigenvector 
#

function B = deflate(A, v)

  [tmp,i] = max(abs(v));
  n = length(A);

  if (i != 1) 
    for k=1:i-1
      for j=1:i-1
	B(k,j) = A(k,j) - (v(k)/v(i))*A(i,j);
      endfor
    endfor
  endif

  if (i!=1 && i!=n)
    for k=i:n-1
      for j=1:i-1
	  B(k,j) = A(k+1,j) - (v(k+1)/v(i))*A(i,j);
	  B(j,k) = A(j,k+1) - (v(j)/v(i))*A(i,k+1);
      endfor
    endfor
  endif

  if (i!=n)
    for k=i:n-1
      for j=i:n-1
	  B(k,j) = A(k+1,j+1) - (v(k+1)/v(i))*A(i,j+1);
      endfor
    endfor
  endif

endfunction
