#!/usr/bin/octave
# Created by Hershal Bhave on 3/01/13
# For M368K HW3, § 9.4 Number 2
# Written in GNU Octave
#
# Description: Finds the symmetric tridiagonal matrix similar to the
# given symmetric matrix
#

function A = householders(A)

  n = length(A);

  for k=1:n-2

    ## Step 2
    q = 0;
    for j=k+1:n
      q+= (A(j,k,k))^2;
    endfor

    ## Step 3
    if A(k+1,k,k) == 0
      alpha = -q^(1/2);
    else 
      alpha = -((q^(1/2)*A(k+1,k,k))/(abs(A(k+1,k,k))));
    endif

    ## Step 4
    rsq = alpha^2 - alpha*A(k+1,k,k);

    ## Step 5
    v(k) = 0; 
    v(k+1) = A(k+1,k,k) - alpha;

    for j=k+2:n
	v(j) = A(j,k,k);
    endfor

    ## Step 6
    for j=k:n
      w(j) = v(j)/sqrt(2*rsq);
    endfor

    P(:,:,k) = eye(size(A,2)) - 2*w'*w;
    A(:,:,k+1) = P(:,:,k)*A(:,:,k)*P(:,:,k);
    
  endfor

endfunction
